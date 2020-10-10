
from os.path import realpath, exists, join
from larch import Group
from larch.fitting import param_group, Parameter
from larch.io import read_xdi
from larch.xafs import feffit_dataset, feffit_transform, feffit, feffit_report, feffpath
from larch.wxlib import _newplot, _plot


def do_fit(self, which):
    if which == 'testrun':
        folder = self.testrun
    else:
        folder = self.baseline
    #end if

    data = read_xdi(join(self.path, 'bromoadamantane.chik'))

    gds = param_group(amp     = Parameter(1.021,     vary=True),
                      enot    = Parameter(4.01,      vary=True),
                      delr    = Parameter(-0.007,    vary=True),
                      brc     = Parameter(expr = '1.9521+delr'),
                      ss      = Parameter(0.003,     vary=True),
                      phir    = Parameter(109.29960 * 3.141592653589793 / 180,   vary=False),
                      cc      = Parameter(1.53780,   vary=False),
                      tanbeta = Parameter(expr = '(brc+cc)*tan(phir/2) / (brc-cc)'),
                      beta    = Parameter(expr = 'atan(tanbeta)'),
                      brc2    = Parameter(expr = '(brc-cc)*cos(phir/2)/cos(beta)'),
                      drh     = Parameter(0.04,      vary=True),
                      ssh     = Parameter(0.005,     vary=True),
                      ss2     = Parameter(expr = 'ss*(brc2/brc)**2'),
                      c3      = Parameter(-0.0007,   vary=True))
    paths = list()
    paths.append(feffpath(realpath(join(folder, "feff0001.dat")),
                          s02    = 'amp',
                          e0     = 'enot',
                          sigma2 = 'ss',
                          deltar = 'delr',
                          third  = 'c3'))
    paths.append(feffpath(realpath(join(folder, "feff0002.dat")),
                          s02    = 'amp',
                          e0     = 'enot',
                          sigma2 = 'ss2',
                          deltar = 'brc2-2.8565'))
    paths.append(feffpath(realpath(join(folder, "feff0003.dat")),
                          s02    = 'amp',
                          e0     = 'enot',
                          sigma2 = 'ssh',
                          deltar = 'drh'))
    paths.append(feffpath(realpath(join(folder, "feff0004.dat")),
                          s02    = 'amp',
                          e0     = 'enot',
                          sigma2 = '(ss+ss2)/2',
                          deltar = '(brc+brc2+cc)/2 - 3.173'))
    paths.append(feffpath(realpath(join(folder, "feff0005.dat")),
                          s02    = 'amp',
                          e0     = 'enot',
                          sigma2 = '(ss+ss2)/2',
                          deltar = '(brc+brc2+cc)/2 - 3.173'))


    trans = feffit_transform(kmin=2.5, kmax=13, kw=(2,1,3), dk=1, window='hanning', rmin=1.25, rmax=3)
    dset  = feffit_dataset(data=data, pathlist=paths, transform=trans)
    fit   = feffit(gds, dset)

    if self.doplot:
        offset = max(dset.data.chir_mag)
        _newplot(dset.data.r,  dset.data.chir_mag+offset, xmax=8, win=2,
                 xlabel=r'$R \rm\,(\AA)$', label='data',
                 ylabel=r'$|\chi(R)| \rm\,(\AA^{-3})$',
                 title='Fit to '+self.folder, show_legend=True)
        _plot(dset.model.r, dset.model.chir_mag+offset, label='fit', win=2)
        _plot(dset.data.r,  dset.data.chir_re, label='data', win=2)
        _plot(dset.model.r, dset.model.chir_re, label='fit', win=2)
    #end if

    if self.verbose:
        print(feffit_report(fit))
    #end if

    return fit
#end def
