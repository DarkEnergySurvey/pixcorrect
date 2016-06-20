#!/usr/bin/env python 
"""
Package of plots to make about sky subtraction.
Called from command line, the first argument is
name of a MiniskyPC file output from sky_pca, and
second argument is the name for a multipage pdf of
diagnostic plots.
"""
import sys, os
import numpy as np
import pylab as pl
from pixcorrect import skyinfo
from argparse import ArgumentParser
cmap = 'cubehelix'

def showMini(mini,vrange=None):
    """
    Make a 2d plot of a mini-sky image.
    Use percentiles of the valid pixels to set the
    color scale if no limit is given.
    """
    if vrange is None:
        vmin = np.percentile(mini.vector(), 1.)
        vmax = np.percentile(mini.vector(), 99.)
    else:
        vmin = -vrange
        vmax = vrange
    pl.imshow(mini.data, interpolation='nearest', origin='lower', aspect='equal',
              vmin=vmin, vmax=vmax, cmap=cmap)
    pl.colorbar()
    return

def showPCA(pcfile):
    """
    Plot all principal components in a row(s)
    """
    pca = skyinfo.MiniskyPC.load(pcfile)
    npc = pca.U.shape[1]
    nrows = (npc-1)/4 + 1
    ncols = min(npc,4)
    fig,axx = pl.subplots(nrows,ncols,squeeze=True)
    fig.set_size_inches(2*ncols,2*nrows)
    for ipc in range(npc):
        irow = ipc/4
        icol = ipc%4
        fig.sca(axx[irow,icol])
        pl.axis('off')
        pc = pca.get_pc(ipc)
        vmin = np.percentile(pc.vector(), 1.)
        vmax = np.percentile(pc.vector(), 99.)
        pl.imshow(pc.data, interpolation='nearest', origin='lower', aspect='equal',
                vmin=vmin, vmax=vmax, cmap=cmap)
        pl.text(10,200,'PC{:d}'.format(ipc),color='white')
    return

def rmsVsPc(pcfile):
    """
    Plot the variance (actually RMS) vs PC number
    """
    rms = skyinfo.MiniskyPC.get_pc_rms(pcfile)
    pl.semilogy(range(len(rms)), rms, 'ro')
    pl.xlim(-0.5,len(rms)+0.5)
    pl.xlabel('PC Number')
    pl.ylabel('RMS signal')
    pl.grid()
    pl.title('RMS vs PC for ' + pcfile)
    return
    
def showResids2(m,model,mask,sh):
    """
    Make an image of residuals in the first 100 frames
    """
    out = np.ones(((sh[0]+2)*10,(sh[1]+2)*10),dtype=float) * -1
    img = np.ones(sh,dtype=float)*-1
    for i in range(10):
        x0=i*(sh[1]+2)
        for j in range(10):
            y0 = j*(sh[0]+2)
            img[mask] = m[:,10*i+j]-model[:,10*i+j]
            print sh, img.shape, x0, y0
            out[y0:y0+sh[0],x0:x0+sh[1]] = img
    return out


def pcaReport(pcafile,pdffile):
    """
    Make a set of plots for quality control on a PCA output
    pcafile is the output of sky_pca
    pdffile is output plots of this routine.
    """
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(pdffile)
    rmsVsPc(pcafile)
    pp.savefig()
    pl.clf()
    showPCA(pcafile)
    pp.savefig()
    pp.close()
    return
