#!/import/www/askap_targets/cgi-bin/python/bin/python
#/import/www/askap_targets/cgi-bin/anaconda2/bin/python
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
#from astroquery import ned
#nedqu = astroquery.ned.Ned.query_region(c,radius=2.*u.arcmin,equinox='J2000.0')
from astropy import coordinates as coord
from astropy import units as u
from astropy.table import Table, Column, vstack
from astropy.io import ascii
from astropy.cosmology import *
#from astropy.cosmology.FlatLambdaCDM import *

#from astropy.coordinates import SkyCoord

from util_html import *
from util_ASKAP import *
from scrape_wise_tiles import *

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
    global NVSSfilename, SUMSSfilename, Srange, nvsslevels, sumsslevels, inS, inN
    nvsslevels=[]
    sumsslevels=[]
    #used for the contours in flux
    
    print "Content-Type: text/html\n\n"
    print_page_head(rootURL + '/style.css', 'Image Display for Radio Sources')
    print "<html><body>"
    print '''
    <div style="width:1200px; border:0px; solid black; text-align: left;
    margin:auto">
    <h2>NVSS & SUMSS Images:</h2>
    '''
    print '''
    <table style="text-align: center; margin:auto;">
    <tbody>
    <tr>
    '''
    print '''
    <td>
    <div class="figureboxborder" style="width:500px">
    '''
    #filename is the name of the image file
    form = cgi.FieldStorage()
    inS = False
    inN = False
    try:
	NVSSfilename = form['Nfile'].value
	Nra_zoom = float(form['NX'].value)
    	Ndec_zoom = float(form['NY'].value)
    	SUMSSfilename = form['Sfile'].value
    	Sra_zoom = float(form['SX'].value)
    	Sdec_zoom = float(form['SY'].value)
    	ra = form['ra'].value
    	dec = form['dec'].value
    	wise = form['wise'].value
	maxNVSS = float(form['N_S'].value)/1000
    	maxSUMSS = float(form['S_S'].value)/1000
    except:
	print '''
        <center><h3>Please check your inputted coordinates or settings - wrong format or illegal characters detected!</h3></center>
        '''
	return
    meanNVSS=maxNVSS/12+0.008#np.median(nvss_data_cube)-0.0001
    meanSUMSS=maxSUMSS/12+0.008
    k_ran = [1.0,1.0,1.4,1.5,1.6,1.8,2.0,2.1,2.2, 2.3, 2.4, 2.4]
    #used for contours - trial and error was involved
    #basically, using the total flux density flux going by the database, and create levels of contours
    #with this as a reference. Works in the vast majority of cases!
    ccc = 0;
    n1 = 6
    n2 = 13.0
    con_string = '</br><em style="font-size:12px">Contour levels (Jy) = %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f</em></br></br>'
    con_nvss = ();
    con_sumss = ();
    if maxNVSS > 10 or maxSUMSS > 10:
        n1 = 11
        meanNVSS=maxNVSS/22 + 0.008
        meanSUMSS=maxSUMSS/22 + 0.008
        n2 = 24.0
        con_string = '</br><em style="font-size:12px">Contour levels (Jy) = %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f</br></br></em>'
        
    for k in range(-1,n1):
        nvsslevels.append(meanNVSS+k*maxNVSS/n2*k_ran[ccc])
        sumsslevels.append(meanSUMSS+k*maxSUMSS/n2*(k_ran[ccc]))
        con_nvss += (nvsslevels[ccc],)
        con_sumss +=(sumsslevels[ccc],)
        ccc+=1;
    #contour range is in Jy/beam.

    #make a figure
    fig = pl.figure(figsize=(10.0, 10.0))    
    label = '%.3f_%.3f' % (float(ra), float(dec))
    
    if NVSSfilename != '' and NVSSfilename != 'null':
        try: 
            #if the file already exists, grab it instead of remaking it
            if os.path.exists('images/NVSS_'+label+'.png'):
                #('NVSS_I/I'+NVSSfilename[1:len(NVSSfilename)]+'.gz'):
                print "<a href='%simages/NVSS_%s.png' target='_blank'>" % (rootURL, label)
                print "<img src='%simages/NVSS_%s.png' height='450' width='450' target='_blank'></a>" % (rootURL, label);
                print '<div class="captionbox" style="text-align: center">'
                if maxNVSS != 0:
                    print '<b style="font-size:12px">Total flux density: %0.2f Jy</b>' % maxNVSS
                    print con_string % con_nvss
                print 'NVSS Image <a'
                #only want the intensity images for NVSS, hence the I + filename(etc) business
                print 'href="%s" ' % ('NVSS_I/I'+NVSSfilename[1:len(NVSSfilename)]+'.gz')
                print '''
                type="image/gz">(FITS FILE)</a>
                </div>
                '''
		inN = True
            else:
                Nurl = ''+'NVSS_I/I'+NVSSfilename[1:len(NVSSfilename)]+'.gz'
                apFig = ap.FITSFigure(Nurl, linewidth=1.0, figure=fig)
                apFig.show_grayscale()#cmap='gist_yarg')#, stretch='arcsinh')
                apFig.set_tick_labels_font(size='small')
                pl.xlim(Nra_zoom - 25, Nra_zoom + 25)
                pl.ylim(Ndec_zoom - 25, Ndec_zoom + 25)
                apFig.show_contour(Nurl, levels=nvsslevels, colors = 'b', alpha=0.5)
                pl.title('NVSS Source')
                #following comments was a temp fix to a problem in saving images
                #Nsio = cStringIO.StringIO()
                #Nurl2 = ''+'images/NVSS_%s.png' % label
                #fig.savefig(Nurl2, format='png')
                #pl.savefig(Nsio, format='png')
                pl.savefig('images/NVSS_'+label+'.png', format='png')
                #print '<img src="data:image/png;base64,%s"/>' % sio.getvalue().encode("base64").strip()
                #msvcrt.setmode(sys.stdout.fileno(), os.O_BINARY) # Needed this on windows, IIS
                #print sys.stdout.write(sio.getvalue())
                print "<a href='%simages/NVSS_%s.png' target='_blank'>" % (rootURL, label)
                print "<img src='%simages/NVSS_%s.png' height='450' width='450' target='_blank'></a>" % (rootURL, label);
                print '<div class="captionbox" style="text-align: center">'
                if maxNVSS != 0:
                    print '<b style="font-size:12px">Total flux density: %0.2f Jy</b>' % maxNVSS
                    print con_string % con_nvss
                print 'NVSS Image <a'
                print 'href="%s" ' % ('NVSS_I/I'+NVSSfilename[1:len(NVSSfilename)]+'.gz')
                print '''
                type="image/gz">(FITS FILE)</a>
                </div>
                
                '''
		inN = True
                
        except Exception:
            print "No such image for NVSS object available."
    else:
        print "This object was not found within NVSS."

    print '''
    </div>
    </td>
    '''
    print '''
    <td>
    <div class="figureboxborder" style="width:500px">
    '''

    #second figure for SUMSS - same structure as before! 
    fig2 = pl.figure(figsize=(10.0, 10.0))
    if str(SUMSSfilename) != '' and str(SUMSSfilename) != 'null':
        #pl.imread('images/test.png') #this works
        try:
            if os.path.exists('images/SUMSS_'+label+'.png'):
                print "<a href='%simages/SUMSS_%s.png' target='_blank'>" % (rootURL, label);
                print "<img src='%simages/SUMSS_%s.png' height='450' width='450' target='_blank'></a>" % (rootURL, label);
                print '<div class="captionbox" style="text-align: center">'
                if maxSUMSS != 0:
                    print '<b style="font-size:12px">Total flux density: %0.2f Jy</b>' % maxSUMSS
                    print con_string % con_sumss
                print 'SUMSS Image <a'
                print 'href="%s" ' % (rootURL+'SUMSS_J/'+SUMSSfilename+'.FITS')
                print '''
                type="image/fits">(FITS FILE)</a>
                </div>
                '''
		inS = True
            else:  
                Surl = ''+'SUMSS_J/'+SUMSSfilename+'.FITS'
                apFig2 = ap.FITSFigure(Surl, linewidth=1.0, figure=fig2)
                apFig2.show_grayscale()#cmap='gray', stretch='arcsinh')
                apFig2.set_tick_labels_font(size='small')
                pl.xlim(Sra_zoom - 25, Sra_zoom + 25)
                pl.ylim(Sdec_zoom - 25, Sdec_zoom + 25)
                apFig2.show_contour(Surl, levels=sumsslevels, colors = 'r', alpha=0.5)
                pl.title('SUMSS Source')
                pl.savefig('images/SUMSS_'+label+'.png', format='png')
                print "<a href='%simages/SUMSS_%s.png' target='_blank'>" % (rootURL, label);
                print "<img src='%simages/SUMSS_%s.png' height='450' width='450' target='_blank'></a>" % (rootURL, label);
                print '<div class="captionbox" style="text-align: center">'
                if maxSUMSS != 0:
                    print '<b style="font-size:12px">Total flux density: %0.2f Jy</b>' % maxSUMSS
                    print con_string % con_sumss
                print 'SUMSS Image <a'
                print 'href="%s" ' % (rootURL+'SUMSS_J/'+SUMSSfilename+'.FITS')
                print '''
                type="image/fits">(FITS FILE)</a>
                </div>
                '''
		inS = True
                         
        except Exception:
            print "No such image for SUMSS object available."

    else:
        print "This object was not found within SUMSS."
    #pl.tight_layout()
    #pl.show();
    
    print '''
    </div>
    </td>
    </tbody>
    </table> 
    '''
    #print '</div>'
    status = 'no'
    #this section is for the wise images, assuming the user wanted them

    if wise == 'true':
        try:
            bands = ['W3.4', 'W4.6', 'W12', 'W22']
            cc = 0;dd = 0;
            #check if all four wise images/fits files already exist - if so, we can skip the process of creating them
            for i in bands:
                im_path = 'images/WISE_%.3f_%.3f_%s.png' % (float(ra), float(dec), i)
                im_path2 = 'images/WISE_%.3f_%.3f_%s.fits' % (float(ra), float(dec), i)
                if os.path.exists(im_path):
                    cc += 1;
                if os.path.exists(im_path2):
                    dd += 1;
            if cc == 4:
                status = 'exists'
            elif dd == 4:
                status = 'no'
                #in this case the images have not been made, but the fits exist, so we don't need to get the fits files
            else:
                status = 'no'
                #here we need to get the fits files first. 
                #getWise is a function in scrape_wise_tiles
                getWise(ra, dec)
        except Exception:
            print '</br></br>Cannot currently get WISE data.'
        plotWise(ra,dec,status)
    #print (datetime.now()-startTime)

    sed_plot(ra,dec)
    print '''
<!-- Start of StatCounter Code for Default Guide -->
<script type="text/javascript">
var sc_project=9930177;
var sc_invisible=1;
var sc_security="b06b31d6";
var scJsHost = (("https:" == document.location.protocol) ?
"https://secure." : "http://www.");
document.write("<sc"+"ript type='text/javascript' src='" + scJsHost+
"statcounter.com/counter/counter.js'></"+"script>");
</script>
<noscript><div class="statcounter"><a title="web analytics"
href="http://statcounter.com/" target="_blank"><img class="statcounter"
src="http://c.statcounter.com/9930177/0/b06b31d6/1/" alt="web
analytics"></a></div></noscript>
<!-- End of StatCounter Code for Default Guide -->
'''
    
    print "</body></html>"


def sed_plot(ra,dec):
    label = '%.3f_%.3f' % (float(ra), float(dec))
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
    <div class="figureboxborder" style="width:660px">
    '''
    if os.path.exists('images/sed_'+label+'.png'):
            print "<a href='%simages/sed_%s.png' target='_blank'>" % (rootURL, label);
            print "<img src='%simages/sed_%s.png' height='450' width='600' target='_blank'></a>" % (rootURL, label);
            print '<div class="captionbox" style="text-align: center">'
            print 'SED from NED'
    else:
            print "<a href='%sgenerate_sed.py?ra=%s&dec=%s' target='_blank'>" % (cgiURL, ra, dec)
            print 'Click to generate the SED from NED!</a>'
            print '</br> A second-order fit will also be attempted for the available points.'
    print '''
    </div>
    </td>
    </tbody>
    </table> 
    '''
    
def plotWise(ra, dec, status='no'):
    #similar to NVSS or SUMSS case
    global NVSSfilename, SUMSSfilename, Srange, nvsslevels, inS, inN
    print '</div>'
    print '''
    <div style="width:1200px; border:0px; solid black; text-align: left;
    margin:auto">
    <h2>WISE Images:</h2>
    '''
    print '''
    <table style="text-align: center; margin:auto;">
    <tbody>
    <tr>
    </tr><tr>
    '''
    bands = ['W3.4', 'W4.6', 'W12', 'W22']
    cc = 0;
    rgbFileLst = [];
    for i in bands:
        fig3 = pl.figure(figsize=(10,10));
        print '''
        <td>
        <div class="figureboxborder" style="width:320px">
        '''        
        try:
            label = '%.3f_%.3f_%s' % (float(ra), float(dec), i)
            if status != 'exists':
                apFig3 = ap.FITSFigure('images/WISE_'+label+'.fits', linewidth=1.0,
                                       figure=fig3)
                #put on NVSS/SUMSS contours if they exist
                if NVSSfilename != '' and NVSSfilename != 'null' and inN != False:
                    apFig3.show_contour('NVSS_I/I'+NVSSfilename[1:len(NVSSfilename)]+'.gz',
                                        levels=nvsslevels, colors = 'b', alpha=1)
                if SUMSSfilename != '' and SUMSSfilename != 'null' and inS != False:
                    apFig3.show_contour('SUMSS_J/'+SUMSSfilename+'.FITS', levels=sumsslevels,
                                        colors='r', alpha=1)
                apFig3.show_colorscale(cmap='gist_yarg', stretch='arcsinh')
                apFig3.set_tick_labels_font(size='small')
                pl.title('WISE ('+i[1:len(i)]+' band)')
                circle_x = (263)/2.
                circle_y = (263)/2.
                #263 is the .fits dimension, so hence want circle at centre
                circ=pl.Circle((circle_x,circle_y), radius=5.5, color='r', fill=False,
                               linewidth=3)
                fig3.gca().add_artist(circ)
               
                #apFig3.show_contour('SUMSS_J/'+SUMSSfilename+'.FITS', levels=Srange,
                #                    colors='r', alpha=0.5)
                pl.savefig('images/WISE_'+label+'.png', format='png')
                           
            print "<a href='%simages/WISE_%s.png' target='_blank'>" % (rootURL, label)
            print "<img src='%simages/WISE_%s.png' width='300'></a>" % (rootURL, label);

            print '''
            <div class="captionbox" style="text-align: center">
            WISE %s microns image <a
            ''' % (i[1:len(i)])
            print 'href="%s" ' % ('WISE_'+label+'.fits')
            print '''
            type="image/fits">(FITS FILE)</a>
            </div>
            '''
        except Exception:
            print "No such WISE image obtainable. "
        print "</div>"
        print "</td>"
        if cc == 1:
            print '''
            </tr></br><tr>
            </tr><tr>
            '''
        if cc < 3:
            rgbFileLst.append('images/WISE_'+label+'.fits');
        cc += 1
    #now we make a colour image, and as before add the contours and etc to it. so coool
    numbr = '%.3f_%.3f_' % (float(ra), float(dec))
    stretchtype = 'arcsinh'
    ap.make_rgb_image(rgbFileLst[::-1], 'images/WISE_'+numbr+'rgb.png', embed_avm_tags=True, stretch_g=stretchtype,stretch_b=stretchtype)
    fig4 = pl.figure(figsize=(10,10));
    
    apFig4 = ap.FITSFigure('images/WISE_'+numbr+'rgb.png', figure=fig4)
    
    apFig4.show_rgb('images/WISE_'+numbr+'rgb.png')
    if NVSSfilename != '' and NVSSfilename != 'null' and inN != False:
        apFig4.show_contour('NVSS_I/I'+NVSSfilename[1:len(NVSSfilename)]+'.gz',
                                        levels=nvsslevels, colors = 'b', alpha=1)
    if SUMSSfilename != '' and SUMSSfilename != 'null' and inS != False:
        apFig4.show_contour('SUMSS_J/'+SUMSSfilename+'.FITS', levels=sumsslevels,
                                        colors='r', alpha=1)
    apFig4.set_tick_labels_font(size='small')
    #pl.title('WISE ('+i[1:len(i)]+' band)')
    circle_x = (263.)/2.
    circle_y = (263.)/2.
    #263 is the .fits dimension, so hence want circle at centre
    circ=pl.Circle((circle_x,circle_y), radius=5.5, color='r', fill=False,
                               linewidth=3)
    fig4.gca().add_artist(circ)
    pl.savefig('images/WISE_'+numbr+'rgb.png')

    print '''
    </tr></br></div></td>
    </tbody>
    </table> 
    '''
    print '''
    <table style="text-align: center; margin:auto;">
    <tbody>
    <tr>
    '''
    print '''</br>
    <td>
    <div class="figureboxborder" style="width:660px">
    '''  
    print "<a href='%simages/WISE_%srgb.png' target='_blank'>" % (rootURL, numbr)
    print "<img src='%simages/WISE_%srgb.png' width='640'></a>" % (rootURL, numbr);

    print '''
            <div class="captionbox" style="text-align: center">
            WISE colour image (blue = 3.4, green = 4.6, red = 12 microns)
            '''
    print ''' </td>
    </div>
    </tbody>
    </table> 
    '''

    print '</div>'

main()
