import numpy as np
from scipy.interpolate import UnivariateSpline as spline
from scipy import stats

class ErrorPropagationSpline(object):
    """
    Does a spline fit, but returns both the spline value and associated uncertainty.
    """
    def __init__(self, x, y, yerr, N=1000, *args, **kwargs):
        """
        See docstring for InterpolatedUnivariateSpline
        """
        yerr=np.clip(yerr,a_min=0.001,a_max=1)
        yy = np.vstack([y + np.random.normal(loc=0, scale=yerr) for i in range(N)]).T
        self._splines = [spline(x, yy[:, i], *args, **kwargs) for i in range(N)]

    def __call__(self, x, *args, **kwargs):
        """
        Get the spline value and uncertainty at point(s) x. args and kwargs are passed to spline.__call__
        :param x:
        :return: a tuple with the mean value at x and the standard deviation
        """
        x = np.atleast_1d(x)
        s = np.vstack([curve(x, *args, **kwargs) for curve in self._splines])
        t=np.arange(10,40,0.01)
        tpeak_ind = np.vstack([np.argmin(curve(t, *args, **kwargs)) for curve in self._splines])
        tpeak=t[np.round(np.median(tpeak_ind))]
        tpeak_err = np.std(np.vstack([t[ind] for ind in tpeak_ind]))
        return (np.median(s, axis=0), np.std(s, axis=0), tpeak_err)

class ErrorPropagationLinear(object):
    #Does a spline fit, but returns both the spline value and associated uncertainty.
    def __init__(self, x, y, yerr, N=1000, *args, **kwargs):
        """
        See docstring for InterpolatedUnivariateSpline
        """
        yerr=np.clip(yerr,a_min=0.001,a_max=1)
        yy = np.vstack([y + np.random.normal(loc=0, scale=yerr) for i in range(N)]).T
        self._splines = [np.polyfit(x, yy[:, i],1) for i in range(N)]

    def __call__(self, x, *args, **kwargs):
        """
        Get the spline value and uncertainty at point(s) x. args and kwargs are passed to spline.__call__
        :param x:
        :return: a tuple with the mean value at x and the standard deviation
        """
        x = np.atleast_1d(x)
        s = np.vstack([np.polyval(curve, x) for curve in self._splines])
        return (np.mean(s, axis=0), np.std(s, axis=0))