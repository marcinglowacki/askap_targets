#!/import/www/askap_targets/cgi-bin/anaconda2/bin/python
#-----------------------------------------------------------------------------#
#usr.bin.python 

import os
import sys
#sys.path.insert(0, '/import/www/askap_targets/cgi-bin/python/lib/python2.7/site-packages/')
#sys.path.append('/usr/lib/python2.7/dist-packages/')
#sys.path.append('/usr/lib/python2.6/site-packages/')
#sys.path.append('/usr/bin/env')
#sys.path.append('/import/www/askap_targets/cgi-bin/anaconda2/pkgs')
#sys.path.append('/usr/local/lib/python2.7/dist-packages/')
#sys.path.append('/import/www/askap_targets/cgi-bin/python/bin/python')

os.environ['HOME'] = '/tmp/'
os.environ['MPLCONFIGDIR']='/tmp/'

from datetime import datetime
fig3startTime = datetime.now()
#I LOVE THE ABOVE LINE 
import cgitb
#import astropy as ast
cgitb.enable()
import string
import re
import urllib
import numpy as np
import matplotlib
#from scipy import *
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
#import StringIO, Image
#import astropy
import aplpy as ap
import sqlite3
import cgi
from cgi import escape
import urllib
import cStringIO
from astroquery import ned
from astropy import coordinates as coord
from astropy import units as u
from astropy.table import Table, Column, vstack
from astropy.io import ascii
from astropy.cosmology import *
#from astropy.cosmology.FlatLambdaCDM import *

#from astropy.coordinates import SkyCoord

from util_html import *
from util_ASKAP import *
#from scrape_wise_tiles import *

# Set up standard LCDM cosmology
#cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
#print FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3).age(0)
cc = 0;
import logging
logging.disable(logging.CRITICAL)
#logging.basicConfig(level=logging.DEBUG)

import warnings

from pylab import *

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

import warnings
warnings.filterwarnings("ignore")

rootURL = open('env_rootURL.txt').read().strip()
cgiURL = open('env_cgiURL.txt').read().strip()

def main():
    global NVSSfilename, SUMSSfilename, Srange, nvsslevels, sumsslevels
    nvsslevels=[]
    sumsslevels=[]
    #used for the contours in flux
    
    print "Content-Type: text/html\n\n"
    print_page_head(rootURL + '/style.css', 'SED creation and fit for Radio Sources')
    print "<html><body>"
    print '''
    <div style="width:1200px; border:0px; solid black; text-align: left;
    margin:auto">
    <h2>Spectral Energy Density:</h2>
    '''
    print '''
    <table style="text-align: center; margin:auto;">
    <tbody>
    <tr>
    '''
    print '''
    <td>
    <div class="figureboxborder" style="width:650px">
    '''

    try:
        form = cgi.FieldStorage()
        ra = form['ra'].value
        dec = form['dec'].value
        label = '%.3f_%.3f' % (float(ra), float(dec))
        c = coord.ICRS(str(ra), str(dec), unit=(u.deg,u.deg))
    except:
        print "Whoops, there's no readable ra and dec fields for this file!</br>"
    try:
        if os.path.exists('images/sed_'+label+'.png'):
            print "<a href='%simages/sed_%s.png' target='_blank'>" % (rootURL, label);
            print "<img src='%simages/sed_%s.png' height='450' width='600' target='_blank'></a>" % (rootURL, label);
            print '<div class="captionbox" style="text-align: center">'
            print 'SED from NED'
        else:
            nedqu = ned.Ned.query_region(c,radius=2.*u.arcmin,equinox='J2000.0')
            dist = np.min(nedqu['Distance (arcmin)'])
            rough_red = ned.Ned.query_object(nedqu['Object Name'][nedqu['Distance (arcmin)']==dist])['Redshift'] 
            try:
                red = ned.Ned.get_table(nedqu['Object Name'][nedqu['Distance (arcmin)']==dist],table='redshifts')
            except:
                pass
            try:
                photo = ned.Ned.get_table(nedqu['Object Name'][nedqu['Distance (arcmin)']==dist],table='photometry')
            except:
                pass
            try:
                z = red['Published Redshift'][0]
        #source['z_ref'] = red['Refcode'][0]
            except:
                z = rough_red
            if np.isnan(float(z)):
                z = -999.
            if (z != -999.):
                pind = 0
                if len(photo['Frequency']) < 2:
                    print 'There is only %d photometrY for this source on NED!\n' % len(photo['Frequency'])
                for j in range(0,len(photo['Frequency'])):
                    if np.isnan(float(photo['NED Photometry Measurement'][pind])):
	       		photo.remove_row(pind)
                    else:
	       		pind += 1
                red = float(z)
                dlum = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3).luminosity_distance(red).value*3.08567758e22
                flux = np.array(photo['NED Photometry Measurement'][:])
                freq = np.array(photo['Frequency'][:])
                freq_prime = freq/(1.+red)
                lum = 4.0*np.pi*np.power(dlum,2)*flux*1.e-26
            # Interpolate flux and luminosity
                flux_spl = np.poly1d(np.polyfit(np.log10(freq[(flux>0.)&(flux<1.e4)&(freq<1.e11)]),
	       					np.log10(flux[(flux>0.)&(flux<1.e4)&(freq<1.e11)]), 2))
                lum_spl = np.poly1d(np.polyfit(np.log10(freq_prime[(lum>0.)&(lum<1.e33)]),
	       					np.log10(lum[(lum>0.)&(lum<1.e33)]), 2))
                lband_freq = np.log10(1.420405752e9/(1.+red))
                lband_flux = np.power(10,flux_spl(lband_freq))
            #source['flux21(Jy)'] = lband_flux
            #print 'Source %s has 21cm flux density = %.3f Jy' \
	    #   		% (source['name'],lband_flux)		
                uv_freq = np.log10(3.29e15)
                uv_lum = np.power(10,lum_spl(uv_freq))
            #print 'Source %s has L_UV = %e W/Hz' \
	    #   		% (source['name'],uv_lum)
            #source['luv(W/HZ)'] = uv_lum
            # Make SED fit plot
                freq_min = 1.e7
                freq_max = 1.e15
                flux_min = 1.e-3
                flux_max = 1.e4
            # lum_min = 1.e10
            # lum_max = 1.e33
                freq_array = np.power(10,np.arange(np.log10(freq_min),np.log10(freq_max),0.1))
                flux_array = np.power(10,flux_spl(np.log10(freq_array)))
            # lum_array = np.power(10,lum_spl(np.log10(freq_array)))
                plt.ioff()
                fig = pl.figure()
                ax = fig.add_subplot(1,1,1)
            # ax.scatter(freq_prime[(lum>0.)&(lum<1.e33)],lum[(lum>0.)&(lum<1.e33)],color='r')
                ax.scatter(freq[(flux>0.)&(flux<1.e4)&(freq<1.e11)],flux[(flux>0.)&(flux<1.e4)&(freq<1.e11)],color='r')		
                ax.scatter(freq[(flux>0.)&(flux<1.e4)&(freq>=1.e11)],flux[(flux>0.)&(flux<1.e4)&(freq>=1.e11)],color='g')	
            # ax.plot(freq_array,lum_array,'--',color='b')
                ax.plot(freq_array,flux_array,'--',color='b')		
            # ax.axvline(x=np.power(10,uv_freq),color=[0.5,0.5,0.5],ls='dashed')
            # ax.axhline(y=1.e23,color=[0.5,0.5,0.5],ls='dashed')
                ax.axvline(x=np.power(10,lband_freq),color=[0.5,0.5,0.5],ls='dashed')
                ax.axhline(y=1.0,color=[0.5,0.5,0.5],ls='dashed')
                ax.text(1.e10,1.e3,r'$z = %.3f$'%(red))
                ax.text(1.e10,3.e2,r'$S_{21} = %.3f\,\mathrm{Jy}$'%(lband_flux))
                ax.text(1.e10,1.e2,r'$L_\mathrm{UV} = %.3e\,\mathrm{W\,Hz^{-1}}$'%(uv_lum))
                pl.suptitle(r'$\mathrm{ra:}$ $\rm{%s,}$ $\rm{dec:}$ $\rm{%s}$'%(c.ra,c.dec))
                ax.set_xscale('log')
                ax.set_yscale('log')
                ax.set_xlabel(r'$\nu\,(\mathrm{Hz})$')
                ax.set_ylabel(r'$S_{\nu}\,(\mathrm{Jy})$')        
            # ax.set_xlabel(r'$\nu\,^{\prime}\,(\mathrm{Hz})$')
            # ax.set_ylabel(r'$L_{\nu\,^{\prime}}\,(\mathrm{W}\,\mathrm{Hz}^{-1})$')        
                ax.set_xlim(freq_min,freq_max)
                ax.set_ylim(flux_min,flux_max)
                pl.savefig('images/sed_'+label+'.png', format='png')

                print "<a href='%simages/sed_%s.png' target='_blank'>" % (rootURL, label);
                print "<img src='%simages/sed_%s.png' height='450' width='600' target='_blank'></a>" % (rootURL, label);
                print '<div class="captionbox" style="text-align: center">'
                print 'SED from NED'
            else:
                print 'This source has no redshift information on NED.\n'
    except:
        print 'This source was either not found within NED, has no redshift information, or an error occurred in reading its photometries.'
    print '''
    </div>
    </td>
    </tbody>
    </table> 
    '''

main()
