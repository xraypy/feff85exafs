import sys

from os.path import realpath, exists, join

from larch import Group
from larch.fitting import param_group, Parameter
from larch.io import read_xdi
from larch.xafs import feffit_dataset, feffit_transform, feffit, feffit_report, feffpath
from larch.wxlib import _newplot, _plot

def do_fit(self, which):

    if which == 'testrun':
        folder = self.testrun
    elif which == 'baseline':
        folder = self.baseline
    else:
        folder = realpath(join(self.folder, 'baseline', which))
    #endif

    print('>>>>>> %s' % folder)

    data = read_xdi(join(self.path, 'Copper.chik'))

    gds = param_group(amp = Parameter(1, vary=True),
                      enot   = Parameter(1e-7,  vary=True),
                      thetad = Parameter(500,   vary=True),
                      temp   = Parameter(10,    vary=False),
                      alpha  = Parameter(1e-7,  vary=True),
                      ss1    = Parameter(0.003, vary=True))


    paths = list()
    for index in range(1, 14):
        nnnn = realpath(join(folder, "feff%4.4d.dat" % index))
        if not exists(nnnn):
            continue
        #end if
        if index > 1:
            sigsqr = 'sigma2_debye(temp, thetad)'
        else:
            sigsqr = 'ss1'
        #end if

        paths.append(feffpath(nnnn, s02='amp', e0='enot',
                              sigma2=sigsqr,
                              deltar='alpha*reff'))

    #end for

    trans = feffit_transform(kmin=3, kmax=16, kw=(2,1,3), dk=1, window='hanning', rmin=1.7, rmax=5.1)
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
    print("Copper Fit Done ", fit)
    return fit
#end def
