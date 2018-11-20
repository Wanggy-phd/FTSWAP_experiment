# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 14:04:53 2018

@author: Jarnd
"""


import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size' : 12})
#rc('text', usetex = True)


def plot_chiandchierror(chi, chi_perror, pnames, title):
    '''
    Plot the real and imaginary part of a chi matrix as a 'city bar' plot.
    The names of the basis in which chi is expressed should be provided in pnames
    pnames should be a 2^(2*n)-long ordered list of the basis names,
    and will be used as x- and yticks.
    Please also provide a title.
    This function can also be used to plot a choi matrix
    '''
    d2 = np.shape(chi)[0]
    _x = np.arange(d2)
    _y = np.arange(d2)
    _xx, _yy = np.meshgrid(_x, _y)
    x, y = _xx.ravel(), _yy.ravel()
    ztot = np.reshape(np.asarray(chi), -1)
    zreal = np.real(ztot)
    zimag = np.imag(ztot)
    bottom = np.zeros_like(zreal)
    width = depth = 0.5
    zlim = [-0.0, np.max(zreal)]
    
    ztoterr = np.reshape(np.asarray(chi_perror), -1)
    zrealerr = np.real(ztoterr)
    zimagerr = np.imag(ztoterr)
    zlimerr = [-0.0, np.max(zrealerr)]


    fig1 = plt.figure(figsize=(12, 11))
    fig1.suptitle(title, fontsize=20)
    ###########################################################################
    ax1 = fig1.add_subplot(211, projection='3d')
    plt.title('Real('+title+')')

    ax1.bar3d(x, y, bottom, width, depth, zreal, shade=True)
    ax1.set_zlim(zlim)

    plt.xticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=55)
    plt.yticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=-60)
#    ax1.set_xlabel(r'$B_{m}$',labelpad=20)
    ax1.set_xlabel(r'$\leftarrow B_{m} \rightarrow$', labelpad=25)
    ax1.set_ylabel(r'$\leftarrow B_{n} \rightarrow$', labelpad=20)
    ax1.xaxis.label.set_fontsize(18)
    ax1.yaxis.label.set_fontsize(18)
    ###########################################################################
    ax2 = fig1.add_subplot(222, projection='3d')
    plt.title('Imag('+title+')')
    ax2.bar3d(x, y, bottom, width, depth, zimag, shade=True)
    ax2.set_zlim(zlim)

    plt.xticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=55)
    plt.yticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=-60)
    ax2.set_xlabel(r'$\leftarrow B_{m} \rightarrow$', labelpad=25)
    ax2.set_ylabel(r'$\leftarrow B_{n} \rightarrow$', labelpad=20)
    ax2.xaxis.label.set_fontsize(18)
    ax2.yaxis.label.set_fontsize(18)
    
    ###########################################################################
    ax2 = fig1.add_subplot(223, projection='3d')
    plt.title('Imag('+title+')')
    ax2.bar3d(x, y, bottom, width, depth, zrealerr, shade=True)
    ax2.set_zlim(zlimerr)

    plt.xticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=55)
    plt.yticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=-60)
    ax2.set_xlabel(r'$\leftarrow B_{m} \rightarrow$', labelpad=25)
    ax2.set_ylabel(r'$\leftarrow B_{n} \rightarrow$', labelpad=20)
    ax2.xaxis.label.set_fontsize(18)
    ax2.yaxis.label.set_fontsize(18)
    ###########################################################################
    ax2 = fig1.add_subplot(224, projection='3d')
    plt.title('Imag('+title+')')
    ax2.bar3d(x, y, bottom, width, depth, zimagerr, shade=True)
    ax2.set_zlim(zlimerr)

    plt.xticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=55)
    plt.yticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=-60)
    ax2.set_xlabel(r'$\leftarrow B_{m} \rightarrow$', labelpad=25)
    ax2.set_ylabel(r'$\leftarrow B_{n} \rightarrow$', labelpad=20)
    ax2.xaxis.label.set_fontsize(18)
    ax2.yaxis.label.set_fontsize(18)
    ###########################################################################
    plt.show()

def plot_city(chi, pnames, title, savename = 'test'):
    '''
    Plot the real and imaginary part of a chi matrix as a 'city bar' plot.
    The names of the basis in which chi is expressed should be provided in pnames
    pnames should be a 2^(2*n)-long ordered list of the basis names,
    and will be used as x- and yticks.
    Please also provide a title.
    This function can also be used to plot a choi matrix
    '''
    d2 = np.shape(chi)[0]
    _x = np.arange(d2)
    _y = np.arange(d2)
    _xx, _yy = np.meshgrid(_x, _y)
    x, y = _xx.ravel(), _yy.ravel()
    ztot = np.reshape(np.asarray(chi), -1)
    zreal = np.real(ztot)
    zimag = np.imag(ztot)
    bottom = np.zeros_like(zreal)
    width = depth = 0.5
    zlim = [-0.0, np.max(zreal)]

    fig1 = plt.figure(figsize=(8, 9))
    fig1.suptitle(title, fontsize=20)
    ###########################################################################
    ax1 = fig1.add_subplot(211, projection='3d')
    plt.title('Real('+title+')')

    ax1.bar3d(x, y, bottom, width, depth, zreal, shade=True)
    ax1.set_zlim(zlim)

    plt.xticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=55)
    plt.yticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=-60)
#    ax1.set_xlabel(r'$B_{m}$',labelpad=20)
    ax1.set_xlabel(r'$\leftarrow B_{m} \rightarrow$', labelpad=15)
    ax1.set_ylabel(r'$\leftarrow B_{n} \rightarrow$', labelpad=15)
    ax1.xaxis.label.set_fontsize(15)
    ax1.yaxis.label.set_fontsize(15)
    ###########################################################################
    ax2 = fig1.add_subplot(212, projection='3d')
    plt.title('Imag('+title+')')
    ax2.bar3d(x, y, bottom, width, depth, zimag, shade=True)
    ax2.set_zlim(zlim)

    plt.xticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=55)
    plt.yticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=-60)
    ax2.set_xlabel(r'$\leftarrow B_{m} \rightarrow$', labelpad=15)
    ax2.set_ylabel(r'$\leftarrow B_{n} \rightarrow$', labelpad=15)
    ax2.xaxis.label.set_fontsize(15)
    ax2.yaxis.label.set_fontsize(15)
    # plt.xlabel('B_{m}',1)
    fig1.subplots_adjust(left = 0.1, right = 0.9,hspace = 0.15, wspace = 0.0, top = 0.95, bottom = 0.15)
    fig1.savefig(fname = 'chiplot-%s.pdf'%(savename), format = 'pdf',bbox_inches = 'tight', pad_inches = 0)
    plt.show()

def plot_city_with_var(chi, chivar, pnames, title):
    '''
    Defunct: bars not stacking but plotting over eachother
    Plot the real and imaginary part of a chi matrix as a 'city bar' plot.
    The standard deviation is plotted as a different-colored bar over the initial chi_matrix.
    So for every m,n there is one bar from chi[m,n]-chivar[m,n],
    and one bar from chi[m,n]-chivar[m,n] to chi[m,n]+chivar[m,n]
    The names of the basis in which chi is expressed should be provided in pnames
    pnames should be a 2^(2*n)-long ordered list of the basis names,
    and will be used as x- and yticks.
    Please also provide a title
    '''
    d2 = np.shape(chi)[0]
    _x = np.arange(d2)
    _y = np.arange(d2)
    _xx, _yy = np.meshgrid(_x, _y)
    x, y = _xx.ravel(), _yy.ravel()
    ztot = np.reshape(np.asarray(chi), -1)
    
        
    zvar = np.abs(np.reshape(np.asarray(chivar), -1))
    zvarreal = np.real(zvar)
    zvarimag = np.imag(zvar)
    
    z1real = abs(np.real(ztot))-zvarreal
    z2real = 2*zvarreal
    z1imag = np.imag(ztot) - zvarimag
    z2imag = 2*zvarimag
    bottom = np.zeros_like(z1real)
    width = depth = 0.5
    zlim = [-0.0, np.max(z1real+z2real)]

    fig1 = plt.figure(figsize=(12, 11))
    fig1.suptitle(title, fontsize=20)
    ###########################################################################
    ax1 = fig1.add_subplot(211, projection='3d')
    plt.title('Real('+title+')')

    ax1.bar3d(x, y, bottom, width, depth, z1real,color=['b'], shade=True, alpha = 1)    
    ax1.bar3d(x,y,z1real, width, depth, z2real,color=['r'], shade = True, alpha = 1)
    
    ax1.set_zlim(zlim)

    plt.xticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=55)
    plt.yticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=-60)
#    ax1.set_xlabel(r'$B_{m}$',labelpad=20)
    ax1.set_xlabel(r'$\leftarrow B_{m} \rightarrow$', labelpad=25)
    ax1.set_ylabel(r'$\leftarrow B_{n} \rightarrow$', labelpad=20)
    ax1.xaxis.label.set_fontsize(18)
    ax1.yaxis.label.set_fontsize(18)
    ###########################################################################
    ax2 = fig1.add_subplot(212, projection='3d')
    plt.title('Imag('+title+')')
    ax2.bar3d(x, y, bottom, width, depth, z1imag, shade=True)
    ax2.set_zlim(zlim)

    plt.xticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=55)
    plt.yticks(np.arange(0.5, d2+0.5, 1), pnames, rotation=-60)
    ax2.set_xlabel(r'$\leftarrow B_{m} \rightarrow$', labelpad=25)
    ax2.set_ylabel(r'$\leftarrow B_{n} \rightarrow$', labelpad=20)
    ax2.xaxis.label.set_fontsize(18)
    ax2.yaxis.label.set_fontsize(18)
    # plt.xlabel('B_{m}',1)

    plt.show()


def plot_contour_fid(pslist,pmlist,height,title,levels=None):
    '''
    Make a contour plot of elements of height on a meshgrid of single and multi qubit error rates.
    The single qubit error rates in pslist are the x-coordinates,
    and the multi qubit error rates in pmlist are the y-coordinate.
    Optional contour levels can be provided in levels.
    Please also provide a title.
    '''
    fig1 = plt.figure(figsize=(12,11))
    fig1.suptitle(title, fontsize=20)
    minlevel = np.around(np.abs(np.min(height)),2)
    maxlevel = np.around(np.abs(np.max(height)),2)
    if levels == None:
        levels = np.linspace(minlevel,maxlevel,num=9)
    CP = plt.contourf(pslist,pmlist,np.real(height),levels,alpha = 0.5)
    Cbar = plt.colorbar(CP)

def plot_surface(pslist,pmlist,height):
    '''
    This function works the same as plot_contour_fid, but makes a surface plot instead.
    
    To do: add option for title
    '''
    fig1 = plt.figure(figsize=(12,11))
    ax1 =fig1.gca(projection='3d')
    X,Y = np.meshgrid(pslist,pmlist)
    surf = ax1.plot_surface(X,Y,height)
    plt.show()
    
def plot_pauli_weight(weights,var, pnames, n):
    '''
    Make a bar plot of the elements in pauli weights in 'weights'. The xticks should be provided in pnames.
    The first bar, corresponding to the weight on the identity matrix, is intentionally left out.
    The error in 'var' is plotted over the normal bars, as in plot_city_with_var
    '''
    fig1 = plt.figure(figsize=(12,11))
    x = np.arange(0,2**(2*n));
    z = weights.copy()
    z1 = z - abs(var)
    z2 = 2*abs(var)
    z1[0] = 0
    z2[0] = 0
    ax1 = fig1.gca()
    ax1.bar(x,z1)
    ax1.bar(x,z2,bottom=z1)
    plt.xticks(np.arange(0,2**(2*n)),pnames)
    