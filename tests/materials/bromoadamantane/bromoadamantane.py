


from os.path import realpath, exists, join

from larch import (Group, Parameter, use_plugin_path)
use_plugin_path('io')
from xdi import read_xdi
use_plugin_path('xafs')
from feffit import feffit_dataset, feffit_transform, feffit, feffit_report
from feffdat import feffpath, _ff2chi
use_plugin_path('wx')
from plotter import (_newplot, _plot)

def do_fit(self, which):

    if which == 'testrun':
        folder = self.testrun
    else:
        folder = self.baseline
    #end if

    data = read_xdi(join(self.path, 'bromoadamantane.chik'), _larch=self._larch)

    gds = Group(amp     = Parameter(1.021,     vary=True,  _larch=self._larch),
                enot    = Parameter(4.01,      vary=True,  _larch=self._larch),
                delr    = Parameter(-0.007,    vary=True,  _larch=self._larch),
                brc     = Parameter(expr = '1.9521+delr',  _larch=self._larch),
                ss      = Parameter(0.003,     vary=True,  _larch=self._larch),
                phir    = Parameter(109.29960 * 3.141592653589793 / 180,   vary=False, _larch=self._larch),
                cc      = Parameter(1.53780,   vary=False, _larch=self._larch),
                tanbeta = Parameter(expr = '(brc+cc)*tan(phir/2) / (brc-cc)', _larch=self._larch),
                beta    = Parameter(expr = 'atan(tanbeta)', _larch=self._larch),
                brc2    = Parameter(expr = '(brc-cc)*cos(phir/2)/cos(beta)', _larch=self._larch),
                drh     = Parameter(0.04,      vary=True,  _larch=self._larch),
                ssh     = Parameter(0.005,     vary=True,  _larch=self._larch),
                ss2     = Parameter(expr = 'ss*(brc2/brc)**2', _larch=self._larch),
                c3      = Parameter(-0.0007,   vary=True,  _larch=self._larch),
                _larch=self._larch  )

    paths = list()
    paths.append(feffpath(realpath(join(folder, "feff0001.dat")),
                          s02    = 'amp',
                          e0     = 'enot',
                          sigma2 = 'ss',
                          deltar = 'delr',
                          third  = 'c3', _larch=self._larch))
    paths.append(feffpath(realpath(join(folder, "feff0002.dat")),
                          s02    = 'amp',
                          e0     = 'enot',
                          sigma2 = 'ss2',
                          deltar = 'brc2-2.8565', _larch=self._larch))
    paths.append(feffpath(realpath(join(folder, "feff0003.dat")),
                          s02    = 'amp',
                          e0     = 'enot',
                          sigma2 = 'ssh',
                          deltar = 'drh', _larch=self._larch))
    paths.append(feffpath(realpath(join(folder, "feff0004.dat")),
                          s02    = 'amp',
                          e0     = 'enot',
                          sigma2 = '(ss+ss2)/2',
                          deltar = '(brc+brc2+cc)/2 - 3.173', _larch=self._larch))
    paths.append(feffpath(realpath(join(folder, "feff0005.dat")),
                          s02    = 'amp',
                          e0     = 'enot',
                          sigma2 = '(ss+ss2)/2',
                          deltar = '(brc+brc2+cc)/2 - 3.173', _larch=self._larch))


    trans = feffit_transform(kmin=2.5, kmax=13, kw=(2,1,3), dk=1, window='hanning', rmin=1.25, rmax=3, _larch=self._larch)
    dset  = feffit_dataset(data=data, pathlist=paths, transform=trans, _larch=self._larch)
    fit   = feffit(gds, dset, _larch=self._larch)

    if self.doplot:
        offset = max(dset.data.chir_mag)
        _newplot(dset.data.r,  dset.data.chir_mag+offset, xmax=8, win=2,
                 xlabel=r'$R \rm\,(\AA)$', label='data',
                 ylabel=r'$|\chi(R)| \rm\,(\AA^{-3})$',
                 title='Fit to '+self.folder, show_legend=True, _larch=self._larch)
        _plot(dset.model.r, dset.model.chir_mag+offset, label='fit', win=2, _larch=self._larch)
        _plot(dset.data.r,  dset.data.chir_re, label='data', win=2, _larch=self._larch)
        _plot(dset.model.r, dset.model.chir_re, label='fit', win=2, _larch=self._larch)
    #end if
    
    if self.verbose:
        print feffit_report(fit, _larch=self._larch)
    #end if

    return fit
#end def
