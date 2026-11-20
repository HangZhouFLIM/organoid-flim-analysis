"""Bi/tri Poisson reconvolution with intact sampled IRF, phase and pulse train.

IRF samples act as fixed quadrature weights at their recorded times. Exponential
emission is integrated over native photon bins, including preceding pulses. The
IRF shape is never estimated from the ROI. Only a common temporal shift is fitted.
"""
import numpy as np
from scipy.optimize import least_squares

def response(amplitudes,taus,background,shift_ns,irf_times_ns,irf_weights,dt_ns,period_ns,n):
    left=np.arange(n)*dt_ns
    width=np.minimum(dt_ns,np.maximum(0,period_ns-left))
    u=(left[:,None]-np.asarray(irf_times_ns)[None,:]-shift_ns)%period_ns
    v=u+width[:,None]
    out=background*width/dt_ns
    for a,t in zip(amplitudes,taus):
        mass=np.where(v<=period_ns,np.exp(-u/t)*(-np.expm1(-width[:,None]/t)),
                      np.exp(-u/t)-np.exp(-period_ns/t)-np.expm1(-np.maximum(v-period_ns,0)/t))
        out=out+a*t*(mass@irf_weights)/(-np.expm1(-period_ns/t))
    return np.maximum(out,1e-12)

def stored_irf(report):
    q=report['irf'];dt=report['bin_width_ns'];n=report['native_bins']
    assert np.isclose(q['sourceStepSeconds']*1e9,dt,rtol=1e-7)
    assert np.isclose(q['referenceDecayStepSeconds']*1e9,dt,rtol=1e-7)
    assert n==int(np.ceil(q['referenceDecayPeriodSeconds']*1e9/dt))
    w=np.asarray(q['sourceValues'],dtype=float)
    if not np.all(np.isfinite(w)) or np.min(w)<0 or w.sum()<=0:
        raise ValueError('Usable matching IRF required')
    origin=q['referenceDecayStartSeconds']-q['referenceDecayStepSeconds']/2
    t=(q['sourceStartSeconds']+np.arange(len(w))*q['sourceStepSeconds']-origin)*1e9
    use=w>0
    return t[use],w[use]/w.sum()

def fit(y,irf_times,irf_weights,dt,period,k,initial_bi=None,*,fit_end=None):
    if k not in (2,3):raise ValueError('Only bi- and tri-exponential fitting is permitted')
    y=np.asarray(y,dtype=float);n=len(y);end=n-1 if fit_end is None else int(fit_end);yy=y[:end]
    if not 2*k+2<end<=n:raise ValueError('Insufficient fitted bins or invalid exclusive end')
    if dt<=0 or period<=0 or not np.isclose(n,np.ceil(period/dt)):raise ValueError('Native grid does not match recorded period')
    assert y.ndim==1 and np.all(y>=0) and np.all(np.isfinite(y))
    assert len(irf_times)==len(irf_weights) and np.isclose(np.sum(irf_weights),1)
    assert yy.sum()>0
    def pred(p):
        return response(np.exp(p[1:k+1]),np.exp(p[k+1:2*k+1]),np.exp(p[0]),p[-1]*dt,irf_times,irf_weights,dt,period,n)
    def residual(p):
        mu=pred(p)[:end];positive=yy>0
        term=mu.copy();u=(mu[positive]-yy[positive])/yy[positive]
        term[positive]=yy[positive]*(u-np.log1p(u))
        return np.sign(mu-yy)*np.sqrt(2*np.maximum(term,0))
    lo=np.log([1e-9,*[1e-9]*k,*[.005]*k]).tolist()+[-3.]
    hi=np.log([y.max(),*[y.max()*1e4]*k,*[20.]*k]).tolist()+[3.]
    guesses=[]
    sets=[(.15,2),(.4,3),(.03,1.7)] if k==2 else [(.08,.6,2.5),(.3,1.2,4),(.1,1,5)]
    for taus in sets:
        for shift in (0.,-1.,-2.):
            guesses.append(np.log([max(1,y[:6].mean()),*[y.max()/dt]*k,*taus]).tolist()+[shift])
    if k==3 and initial_bi:
        a=initial_bi['amplitudes'];t=initial_bi['lifetimes_ns']
        guesses.append(np.log([max(1e-8,initial_bi['baseline_counts_per_bin']),a[0],.8*a[1],.2*a[1],t[0],t[1],min(15,2*t[1])]).tolist()+[initial_bi['irf_shift_bins']])
    trials=[];best=None
    for p in guesses:
        sol=least_squares(residual,np.clip(p,np.array(lo)+1e-8,np.array(hi)-1e-8),bounds=(lo,hi),
                          max_nfev=2000,ftol=1e-10,xtol=1e-10,gtol=1e-7)
        trials.append(dict(success=bool(sol.success),status=int(sol.status),nfev=int(sol.nfev),deviance=float(2*sol.cost)))
        if np.all(np.isfinite(sol.x)) and np.isfinite(sol.cost) and (best is None or sol.cost<best.cost):best=sol
    if best is None:raise RuntimeError('All fitting initializations failed')
    p=best.x;prediction=pred(p);a=np.exp(p[1:k+1]);t=np.exp(p[k+1:2*k+1]);order=np.argsort(t);a=a[order];t=t[order]
    pearson=(yy-prediction[:end])/np.sqrt(prediction[:end]);chi=float(pearson@pearson);dof=end-(2*k+2)
    lag=float(np.corrcoef(pearson[:-1],pearson[1:])[0,1]);fractions=a/a.sum();flags=[]
    if not best.success:flags.append('optimizer_not_converged')
    if chi/dof>2:flags.append('high_reduced_chisq_diagnostic_threshold_2')
    if abs(lag)>.3:flags.append('structured_residuals_lag1_abs_gt_0_3')
    if min(t)<dt:flags.append('component_shorter_than_native_bin')
    if np.any(t[1:]/t[:-1]<1.2) or min(fractions)<.01:flags.append('weakly_separated_or_small_component')
    if min(t)<.0051 or max(t)>19.99:flags.append('lifetime_at_search_boundary')
    if abs(p[-1])>2.99:flags.append('irf_shift_at_search_boundary')
    if np.exp(p[0])<1e-4:flags.append('background_at_zero_boundary')
    return dict(status='finite_fit',model_components=k,fit_start_bin_zero_based=0,fit_end_bin_exclusive=end,
                degrees_of_freedom=dof,raw_chisq=chi,reduced_chisq=chi/dof,poisson_deviance=float(2*best.cost),
                baseline_counts_per_bin=float(np.exp(p[0])),amplitudes=a.tolist(),amplitude_units='counts per ns at excitation',
                lifetimes_ns=t.tolist(),amplitude_fractions=fractions.tolist(),intensity_fractions=(a*t/(a@t)).tolist(),
                amplitude_weighted_lifetime_ns=float(a@t/a.sum()),intensity_weighted_lifetime_ns=float(a@(t*t)/(a@t)),
                irf_shift_bins=float(p[-1]),irf_shift_ns=float(p[-1]*dt),residual_lag1_correlation=lag,
                converged=bool(best.success),optimizer_message=best.message,optimizer_optimality=float(best.optimality),
                finite_initialisations=len(trials),trials=trials,jacobian_singular_values=np.linalg.svd(best.jac,compute_uv=False).tolist(),
                qc_flags=flags,qc_pass=not flags,
                interpretation='Technical goodness-of-fit only; inspect IRF provenance; calculated IRFs are not independent calibrations and components have no molecular assignment'),prediction
