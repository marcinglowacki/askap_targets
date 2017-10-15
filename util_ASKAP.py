#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     Previously util_CORNISH.py                                        #
#                                                                             #
# PURPOSE:  Common functions first used for the CORNISH CGI server scripts.   #
#                                                                             #
# MODIFIED: 06-May-2013 by C. Purcell, changed by M. Glowacki (mid 2014)      #
#                                                                             #
# TODO:                                                                       #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# get_base_cat_sql      ... return SQL to get the published catalogue         #
# query_catalogue_cone  ... perform a cone search on the catalogue            #
# query_catalogue_name  ... query the catalogue based on a name regexp        #
# query_catalogue_radec ... query the catalogue inside a RA-Dec box           #
# query_pnt_in_tiles    ... query the tileIDs containing a point              #
# query_tiles_with_pnt  ... query the tileIDs containing a point              #
# query_uvdata_with_pnt ... query the split uvdata contaiing a point          #
# print_cat_html_alt    ... print HTML cat with options to reload as ascii    #
# print_cat_tab_html    ... print the HTML catalogue table                    #
# print_cat_tab_ascii   ... print the ASCII catalogue table                   #
# print_cat_tab_csv     ... print the CSV catalogue table                     #
# print_cat_tab_kvis    ... print the KVIS catalogue table                    #
# print_cat_tab_aips    ... print the AIPS catalogye table                    #
# print_cat_tab_ds9     ... print the DS9 catalogue table                     #
# parse_coord_string    ... parse a user-input coordinate string              #
# parse_coord_file_cone ... parse the user-input cone-search catalogue        #
# parse_coord_file_cutter . pare the user-input cutout catalogue              #
# parse_radius          ... parse the user-input radius/size string           #
# get_clicked_coord     ... return world-coordinate clicked on simple image   #
# check_type            ... identify the query type for .CSV data output      #
#=============================================================================#

#note - most of these are not used, left over from the CORNISH work
#used: get_base_cat_sql, query_catalogue_cone, and printing for html/csv.
#these functions have been changed accordingly for my database. 

import os
import sys
import re
import math as m

#import Image as im
#import pyfits as pf
#import pywcs as pw
#from pyslalib import slalib
import sqlite3

from util_html import *
#from util_misc import *
#from util_fits import *

#conn = sqlite3.connect('test_db.db')
#cursor = conn.cursor()

#-----------------------------------------------------------------------------#
# Basic catalogue SQL query to retreive the published catalogue               #
#-----------------------------------------------------------------------------#
def get_base_cat_sql():
    #I changed this, gets the more important columns from the base NVSS/SUMSS side, with AT20G data too
    selStat = '''SELECT  
    t1.N_St,
    t1.N_e_St,
    t1.CENTER_X,
    t1.CENTER_Y,
    t1.FIELD,
    t1.tot_flux_dens,
    t1.tot_flux_dens_uncer,
    t1.mosaic_name, 
    t1.mosaic_x,  
    t1.mosaic_y,
    t1.AT20G_Name,
    t1.S_20Ghz,
    t1.S_20Ghz_e,
    t1.S_8Ghz,
    t1.S_8Ghz_e,
    t1.S_5Ghz,
    t1.S_5Ghz_e 
    '''
#    selStat = '''SELECT *'''
    fromStat = '''FROM complete_all t1'''
    return selStat, fromStat


#-----------------------------------------------------------------------------#
# Perform a cone-search query on the database                                 #
#-----------------------------------------------------------------------------#
def query_catalogue_cone(ra_deg, dec_deg, rad_deg):
    #changed as sqlite does not have maths functions, or at least could not find them at the time
    #commented out MySQL code, used pythagoras theorem instead
    selSQL, fromSQL = get_base_cat_sql()
    ## selSQL = """,  
    ## degrees(acos(cos(radians(90-%s))
    ## *cos(radians(90-dec_new))
    ## +sin(radians(90-%s))
    ## *sin(radians(90-dec_new))
    ## *cos(radians(ra_new-%s))))
    ## AS dist
    ## """ % (dec_deg,
    ##        dec_deg,     
    ##        ra_deg)
    posList = (float(ra_deg), float(ra_deg), float(dec_deg), float(dec_deg))
    selSQL2 = ', ((ra_new - %f)*(ra_new - %f) + (dec_new - %f)*(dec_new - %f)) AS dist ' % posList
    
    sql = selSQL + ', ra_new, dec_new ' + selSQL2 + fromSQL
    sql += ' WHERE dist < (%f*%f)' % (float(rad_deg), float(rad_deg))
    ## sql += """ 
    ## WHERE
    ## degrees(acos(cos(radians(90-%s))
    ## *cos(radians(90-deg_new))
    ## +sin(radians(90-%s))
    ## *sin(radians(90-deg_new))
    ## *cos(radians(ra_new-%s)))) < %s
    ## """ % (dec_deg, dec_deg, ra_deg, rad_deg)
    ##sql += "ORDER by dist"
    sql = sql.replace('\n', ' ')

    #rowLst = select_into_namedlst(cursor, sql)

    return sql
   
#-----------------------------------------------------------------------------#
# Query the sources matching a name string or regular expression              #
#-----------------------------------------------------------------------------#
def query_catalogue_name(cursor, srchStr, sigmaLimit=7.0, retArti=True):
    
    selSQL, fromSQL = get_base_cat_sql()
    sql = selSQL + fromSQL
    sql += """ 
    WHERE
    catalogue_master.sourceName REGEXP %s
    AND catalogue_master.photSigma >= %s
    """
    if not retArti:
        sql += '''
        AND catalogue_id_master.artifactProb = "Unlikely"
        '''
    rowLst = select_into_namedlst(cursor, sql, [srchStr, sigmaLimit])

    return rowLst


#-----------------------------------------------------------------------------#
# Query the source notes matching a string or regular expression              #
#-----------------------------------------------------------------------------#
def query_catalogue_desc(cursor, srchStr, sigmaLimit=7.0, retArti=True):
    
    selSQL, fromSQL = get_base_cat_sql()
    fromSQL += """
    LEFT JOIN catalogue_id_notes
    ON catalogue_master.sourceName = catalogue_id_notes.sourceName
    """
    
    sql = selSQL + fromSQL
    sql += """ 
    WHERE
    catalogue_id_notes.idNotes REGEXP %s
    AND catalogue_master.photSigma >= %s
    """
    if not retArti:
        sql += '''
        AND catalogue_id_master.artifactProb = "Unlikely"
        '''
    rowLst = select_into_namedlst(cursor, sql, [srchStr, sigmaLimit])

    return rowLst



#-----------------------------------------------------------------------------#
# Query the sources within a RA & Dec bounded region                          #
#-----------------------------------------------------------------------------#
def query_catalogue_radec(cursor, raMax_deg, raMin_deg, decMax_deg, decMin_deg,
                          sigmaLimit=7.0):
    
    selSQL, fromSQL = get_base_cat_sql()
    sql = selSQL + fromSQL
    sql += """ 
    WHERE catalogue_master.photSigma >= %s
    AND catalogue_master.RAwcent_hr*15 >= %s
    AND catalogue_master.RAwcent_hr*15 <= %s
    AND catalogue_master.DecWcent_deg >= %s
    AND catalogue_master.DecWcent_deg <= %s
    """
    rowLst = select_into_namedlst(cursor, sql, [sigmaLimit, raMin_deg,
                                                raMax_deg, decMin_deg,
                                                decMax_deg])
    return rowLst

#-----------------------------------------------------------------------------#
# Query the tiles containing a point                                          #
#-----------------------------------------------------------------------------#
def query_pnt_in_tiles(cursor, ra_deg, dec_deg, size_deg):

    tileIDList = []
    radius_deg = abs(size_deg/2.0)
    
    # Query the database for all tiles containing the point clicked
    # Results (ID,RA,DEC,Offset) are ordered by offset.
    query="""
    SELECT DISTINCT tile_id,
    degrees(acos(cos(radians(90.0-%s))
    *cos(radians(90.0-DEC_deg))
    +sin(radians(90.0-%s))
    *sin(radians(90.0-DEC_deg))
    *cos(radians(RA_hr*15.0-%s))))
    AS offset
    FROM tile_coords
    WHERE DEC_deg  <= %s + ((npix_y*pixscale_y/(2.0*3600.0))+%s)
    AND   DEC_deg  >= %s - ((npix_y*pixscale_y/(2.0*3600.0))+%s)
    AND RA_hr*15.0 <= %s + ((npix_x*pixscale_x/(2.0*3600.0))+%s)
    /cos(radians(%s))
    AND RA_hr*15.0 >= %s - ((npix_x*pixscale_x/(2.0*3600.0))+%s)
    /cos(radians(%s))
    ORDER BY offset DESC""" % \
    (dec_deg, dec_deg, ra_deg,
     dec_deg, radius_deg,
     dec_deg, radius_deg,
     ra_deg, radius_deg, dec_deg,
     ra_deg, radius_deg, dec_deg)
    cursor.execute(query)
    
    rows = cursor.fetchall()
    if rows:
        for row in rows:
            tileIDList.append(int(row[0]))
    return tileIDList


#-----------------------------------------------------------------------------#
# Query the tiles containing a point                                          #
#-----------------------------------------------------------------------------#
def query_tiles_with_pnt(cursor, ra_deg, dec_deg):

    tileIDList = []
    ra_hr = ra_deg / 15.0
    
    # Query the database for all tiles containing the point clicked
    query="""
    SELECT DISTINCT
    tile_id
    FROM tile_coords
    WHERE DEC_deg <= (%s+((npix_y/2)*(pixscale_y/3600)))
    AND DEC_deg >= (%s-((npix_y/2)*(pixscale_y/3600)))
    AND RA_hr <= (%s+((npix_y/2)*(pixscale_y/(3600*15))))
    AND RA_hr >= (%s-((npix_y/2)*(pixscale_y/(3600*15))))
    """ % (dec_deg, dec_deg, ra_hr, ra_hr)
    cursor.execute(query)
    rows = cursor.fetchall()
    if rows:
        for row in rows:
            tileIDList.append(int(row[0]))
    return tileIDList


#-----------------------------------------------------------------------------#
# Query the tiles containing a point                                          #
#-----------------------------------------------------------------------------#
def query_uvdata_with_pnt(cursor, ra_deg, dec_deg):
    
    # Query the database for all fields containing the point clicked
    sql = """
    SELECT DISTINCT
    observations.id,
    observations.field_name,
    observations.duration,
    observations.sample_time,
    field_coords.RA_hr,
    field_coords.DEC_deg,
    field_coords.l_deg,
    field_coords.b_deg
    FROM observations INNER JOIN field_coords
    ON observations.field_name = field_coords.field_name
    WHERE observations.calcode = ''
    AND observations.duration IS NOT NULL
    AND observations.field_name NOT REGEXP 'B'
    AND observations.useflag = 1
    AND degrees(acos(cos(radians(90-%s))
    *cos(radians(90-field_coords.DEC_deg))
    +sin(radians(90-%s))
    *sin(radians(90-field_coords.DEC_deg))
    *cos(radians(field_coords.RA_hr*15-%s)))) < 7.4/60
    ORDER BY BINARY observations.field_name
    """
    
    rowLst = select_into_namedlst(cursor, sql, [dec_deg, dec_deg, ra_deg])
    
    return rowLst


#-----------------------------------------------------------------------------#
# Generate the HTML code to print the catalogue table                         #
#-----------------------------------------------------------------------------#
def print_cat_html_alt(rootURL, cgiURL, rowLst, showSep=False, showErrs=False, showCal=False, showCom=False, showBers=False, showWP=False):
    
    # Create link to other versions of the table
    scriptName =  os.path.split(os.environ['SCRIPT_FILENAME'])[-1]
    htmlQuery = os.environ['QUERY_STRING']
    asciiQuery = htmlQuery.replace('HTML', 'ASCII')
    asciiQuery = cgiURL + "/" + scriptName + "?" + asciiQuery
    csvQuery = htmlQuery.replace('HTML', 'CSV')
    csvQuery = cgiURL + "/" + scriptName + "?" + csvQuery
    kvisQuery = htmlQuery.replace('HTML', 'KVIS')
    kvisQuery = cgiURL  + "/" + scriptName + "?" + kvisQuery
    ds9Query = htmlQuery.replace('HTML', 'DS9')
    ds9Query = cgiURL  + "/" + scriptName + "?" + ds9Query
    aipsQuery = htmlQuery.replace('HTML', 'AIPS')
    aipsQuery = cgiURL  + "/" + scriptName + "?" + aipsQuery
    print '''
    <h3>Alternative Formats:</h3>
    <p>Click on the following links to reload this page in an 
    alternative format.</br>
    '''
    ## print '<a href="%s">ASCII Format</a> | ' % asciiQuery
    print '<a href="%s">CSV Format</a> ' % csvQuery
    ## print '<a href="%s/downloads/CORNISH_cat_format.txt">(Format Notes)</a><br>' \
    ##       % rootURL
    ## print 'Annotation Files: <a href="%s">KVIS</a>' % kvisQuery
    ## print ' | <a href="%s">DS9</a>' % ds9Query
    ## print ' | <a href="%s">AIPS</a></p></br>' % aipsQuery
    
    # Print the catalogue table
    print '''
    <p style="text-align:center;">Click on a source number to view the
    NVSS and SUMSS images for each source.</br></p>
    <p>If you wish to view the WISE images as well, click the link 'With WISE' for the
    desired source.</br>
    WISE image data obtained from
    <a href="http://hachi.ipac.caltech.edu:8080/montage/">IPAC's image mosaic service</a>, 
    and will <b>add at least a minute to the script run time</b> before
    <em>any</em> images are displayed.</br>
    <strong>The IPAC image mosaic service appears to be undergoing issues, so WISE images may not be able to be currently created.</strong>
Further data to come as well! </br></p>
    
    '''
    #showCom = False;
    #showCom = True;
    
    print_cat_tab_html(cgiURL, rowLst, showSep, showErrs, showCal, showCom, showBers, showWP)

    
#-----------------------------------------------------------------------------#
# Generate the HTML code to print the catalogue table                         #
#-----------------------------------------------------------------------------#
def print_cat_tab_html(cgiURL, rowLst, showSep=False, showErrs=False, showCal=False, showCom=False, showBers=False, showWP=False):
       
    # Print the table header
    showNed = True;
    print '''
    <table class="fancy" style="margin-left:auto; margin-right:auto;
    text-align:center;">
    <tbody>
    <tr>
    <th>No.</th>
    <th>View with</br>WISE data</th>
    <th>RA </br>(J2000)</th>
    <th>Dec </br>(J2000)</th>
    <th>NVSS Flux</br>(mJy)</th>
    <th>NVSS Flux</br>Error (mJy)</th>
    <th>SUMSS Flux </br>Density (mJy)</th>
    <th>SUMSS Flux </br>Error (mJy)</th>
    <th>AT20G FLUX</br> 20 GHz (mJy)</th>
    <th>AT20G Flux 20 Ghz</br>Error (mJy)</th>
    '''
    raExt = 'ra_new'
    decExt = 'dec_new'
    
    if showSep:
        print '<th>Ang. Sep.</br>(arcsec)</th>'
    if showCal:
        print '<th>ATCA Calibrator</br> Name</th>'
        raExt = 'ra_new_Cal'
        decExt = 'dec_new_Cal'
    if showCom:
        print '<th>Alpha</th>'
        print '<th>Reduced </br>Chi-square</th>'
        print '<th>SED Fit </br>PDF File</th>'
        raExt = 'ra_new_Com'
        decExt = 'dec_new_Com'
    if showBers:
        print '<th>BERS source</br>Name</th>'
        print '<th>3C Source</br>Name</th>'
        print '<th>Source</br>Class</th>'
        print '<th>Alpha</th>'
        print '<th>Redshift</th>'
        #unnecessary now that I've done the 'mixed' query separately
        if not showCal:
            raExt = 'ra_new_B'
            decExt = 'dec_new_B'
    if showWP:
        print '<th>WP85 source</br>Name</th>'
        print '<th>Alt. source</br>Name</th>'
        print '<th>WP85 Flux</br>(Jy)</th>'
        print '<th>Class</th>'
        print '<th>Alpha</th>'
        print '<th>Redshift</th>'
        raExt = 'ra_new_WP'
        decExt = 'dec_new_WP'
    if showNed:
        print '<th>Search</br>NED</th>'
        #this is always true, but in case I do not want NED links in future...
    if not showSep:
        print '<th>Within</br>5 degrees</th>'

    print '</tr>'
    
    # Print the source entries
    i = 0
    for e in rowLst:
        i+=1

        NVSSflux = e['N_St']
        NF_hold = NVSSflux
        if NVSSflux == '':
            NF_hold = 0
        NVSSfluxerr = e['N_e_St']
        SUMSSflux = e['tot_flux_dens']
        SF_hold = SUMSSflux
        if SUMSSflux == '':
            SF_hold = 0
        SUMSSfluxerr = e['tot_flux_dens_uncer']
        AT_20Ghz = e['S_20Ghz']
        AT_20Ghzerr = e['S_20Ghz_e']
        if e['FIELD'] == '':
            e['FIELD'] = 'null'
            e['CENTER_X'] = 0
            e['CENTER_Y'] = 0
        if e['mosaic_name'] == '':
            e['mosaic_name'] = 'null'
            e['mosaic_x'] = 0
            e['mosaic_y'] = 0

        #link = cgiURL + 'image_display.py?ra=' + \
        #       str2cgi("%f" % e['ra_new']) + str2cgi("&") + 'dec='+ \
        #       str2cgi("%f" % e['dec_new'])
        #THIS IS A PLACEHOLDER in attempt for second sql search
        #in the image display file, but let's skip that:
        #was testing stuff earlier
        try:
            extra = '?Nfile=' + \
          str2cgi("%s" % e['FIELD']) + '&NX=' + \
          str2cgi("%f" % e['CENTER_X']) +'&NY=' + \
          str2cgi("%f" % e['CENTER_Y']) + '&Sfile=' + \
          str2cgi("%s" % e['mosaic_name']) + '&SX=' + \
          str2cgi("%f" % e['mosaic_x']) + '&SY=' + \
          str2cgi("%f" % e['mosaic_y']) + '&ra=' + \
          str2cgi("%f" % e[raExt]) + '&dec=' + \
          str2cgi("%f" % e[decExt]) + '&N_S=' + \
          str2cgi("%s" % NF_hold) + '&S_S=' + \
          str2cgi("%s" % SF_hold);
        except Exception:
            extra = 'WHOOPS'
        link = cgiURL + 'image_display.py' + extra
        print '<tr>'
        print '<td><a href="%s" target="_blank">%d</a></td>' % (link+'&wise=false', (i))
        print '<td><a href="%s" target="_blank">%s</a></td>' % (link+'&wise=true', ('With WISE'))
        print '<td>%s</td>' % deg2hms(e[raExt])
        print '<td>%s</td>' % deg2dms(e[decExt], doSign=True)
        try:
                print '<td>%.2f</td>' % NVSSflux
                print '<td>%.2f</td>' % NVSSfluxerr
        except Exception:
                print '<td></td>'
                print '<td></td>'
        try:
                print '<td>%.2f</td>' % SUMSSflux
                print '<td>%.2f</td>' % SUMSSfluxerr
        except Exception:
                print '<td></td>'
                print '<td></td>'
        try:
                print '<td>%.2f</td>' % AT_20Ghz
        except Exception:
                print '<td></td>'
        try:
                print '<td>%.2f</td>' % AT_20Ghzerr
        except Exception:
                print '<td></td>'
                
        if showSep:
                try:
                        print '<td>%.2f</td>' % (m.sqrt(float(e['dist'])) * 3600.0)
                except Exception:
                        print '<td>N/A</td>'

        if showCal:
                calLink = 'http://www.narrabri.atnf.csiro.au/calibrators/'
                #calLink += 'calibrator_database_viewcal?source=%s' % e['C_Name']
                try:
                        calLink += 'calibrator_database_viewcal?source=%s' % e['C_Name']
                        print '<td><a href="%s" target="_blank">%s</a></td>' % (calLink, e['C_Name'])
                except Exception:
                        print '<td>N/A</td>'

        if showCom:
                comLink = 'http://www.physics.usyd.edu.au/askap_targets/cgi-bin/images/figures/'
                try:
                        print '<td>%.2f</td>' % (e['alpha'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        print '<td>%.2f</td>' % (e['redchisq'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        comLink += '%s.pdf' % e['NVSS_Name'] 
                        print '<td><a href="%s" target="_blank">%s</a></td>' % (comLink, e['NVSS_Name'])
                except Exception:
                        print '<td>N/A</td>'

        if showBers:
                redshift = str(e['z']);
                if redshift == '0.0':
                        redshift = ''
                bersLink = cgiURL + 'bers_query.py?title=BERSFlux&source=%s' % (e['Source_Name_1Jy'])
                try:
                        print '<td><a href="%s" target="_blank">%s</a></td>' % (bersLink, e['Source_Name_1Jy'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        print '<td>%s</td>' % (e['ThreeC_Name'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        print '<td>%s</td>' % (e['Class'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        print '<td>%.2f</td>' % (e['alpha'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        print '<td>%s</td>' % redshift
                except Exception:
                        print '<td>N/A</td>'

        if showWP:
                redshift = str(e['z']);
                if redshift == '0.0':
                        redshift = ''
                #bersLink = cgiURL + 'bers_query.py?title=BERSFlux&source=%s' % (e['Source_Name_1Jy'])
                try:
                        print '<td>%s</td>' % (e['IAUname'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        print '<td>%s</td>' % (e['AltName'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        print '<td>%s</td>' % (e['S11cm'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        print '<td>%s</td>' % (e['Otype'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        print '<td>%.2f</td>' % (e['Sp_Index'])
                except Exception:
                        print '<td>N/A</td>'
                try:
                        print '<td>%s</td>' % redshift
                except Exception:
                        print '<td>N/A</td>'

        if showNed:
                try:
                        nedLink = 'http://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys=Equatorial'
                        nedLink += '&in_equinox=J2000.0'
                        nedLink += '&lon=%s&lat=%s' % (deg2hms(e[raExt]), deg2dms(e[decExt],doSign=True))
                        nedLink += '&radius=2.0&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&search_type='
                        nedLink += 'Near+Position+Search&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z'
                        nedLink += '&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0'
                        nedLink += '&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0'
                        nedLink += '&list_limit=5&img_stamp=YES'

                        print '<td><a href="%s" target="_blank">Search NED</a></td>' % nedLink
                except Exception:
                        print '<td>N/A</td>'
                searchLink = cgiURL + 'cat_query.py?ra=%s&dec=%s&rad=300&title=Pos_Search' % (e[raExt],e[decExt])
                if not showSep:
                    print '<td><a href="%s" target="_blanks">Search</td>' % searchLink
                        
        print '</tr>'
            
    # Closing tags
    print '''
    </tbody>
    </table>
    </br></br>
    '''


#-----------------------------------------------------------------------------#
# Serve the catalogue query results in plain text format                      #
#-----------------------------------------------------------------------------#
def print_cat_tab_ascii_old(rowLst, showSep=False):
                            
    # Print the column headers                      
    print "NAME" + ' '*13,
    print "l" + " "*8,
    print "b" + " "*7,
    print "RA" + " "*9,
    print "DEC" + " "*9,
    print "ANG" + " "*2,
    print "PEAK" + " "*5,
    print "INTEG" + " "*4,
    print "RMS" + " "*6,
    print "AVG" + " "*5,
    print "SIGMA" + " "*1,
    print "CODE",
    if not showSep:
        print '\n',
    else:
        print 'ANGSEP' + " "*3,
        print "| "
    
    print " "*17,
    print "(Deg)" + " "*4,
    print "(Deg)" + " "*3,
    print "(J2000)" + " "*4,
    print "(J2000)" + " "*5,
    print "(\")" + " "*2,
    print "(mJy/bm)" + " "*1,
    print "(mJy)" + " "*4,
    print "(mJy/bm)" + " "*1,
    print "(mJy/bm)" + " "*12,
    if not showSep:
        print '\n',
    else:
        print '(")' + " "*6,
        print "|"

    # Print the catalogue entries
    for e in rowLst:
        print e['sourceName'],
        print '%s' % num2str(e['l_deg'],3,4,0,0,1),
        print '%s ' % num2str(e['b_deg'],3,4,0,0,1),
        print '%s' % deg2dms(e['RAwcent_hr']), 
        print '%s' % deg2dms(e['DecWcent_deg'], doSign=True),
        print '%s' % num2str(e['mergeAngscaleAbs_asec'],2,2,0,0,1),
        print '%s' % num2str(e['mergePeakIntAbs_mJybm'],5,3,0,0,1),
        print '%s' % num2str(e['mergeFluxAbs_mJy'],5,3,0,0,1),
        print '%s' % num2str(e['photMADFMsky_mJybm'],4,3,0,0,1),
        print '%s' % num2str(e['photMedianSky_mJybm'],4,3,0,0,1),
        print '%s ' % num2str(e['photSigma'],4,1,0,0,1),
        print '%s' % e['mType'],
        print ' '*2,
        if showSep:
            print '%s ' % num2str(e['dist'] * 3600.0,6,1,0,0,1),
            print "| ",
        if e.has_key('inCatDict'):
            print '%s ' % num2str(e['inCatDict']['indx'],4,0,0,0,1),
            print ' '.join(e['inCatDict']['extra'])
        else:
            print '\n',


#-----------------------------------------------------------------------------#
# Serve the catalogue query results in plain text format                      #
#-----------------------------------------------------------------------------#
def print_cat_tab_ascii(rowLst, showSep=False, dlz=0):
                            
    unresolved_asec = 1.8
    
    # Print the column headers     
    line = ''
    line += '#Name' + ' '*13
    line += 'l_deg' + ' '*5
    line += 'b_deg' + ' '*5
    line += 'RA_deg' + ' '*4
    line += 'Dec_deg' + ' '*3
    line += 'dRA_asec' + ' '*1
    line += 'dDec_asec' + ' '*1
    line += 'Peak_mJybm' + ' '*1
    line += 'dPeak_mJybm' + ' '*1
    line += 'Flux_mJy' + ' '*1
    line += 'dFlux_mJy' + ' '*1
    line += 'Angscale_asec' + ' '*1
    line += 'dAngscale_asec' + ' '*1
    line += 'AngscaleDecon_asec' + ' '*1
    line += 'gaussMajor_asec' + ' '*1
    line += 'dGaussMajor_asec' + ' '*1
    line += 'gaussMinor_asec' + ' '*1
    line += 'dGaussMinor_asec' + ' '*1
    line += 'gaussPosangle_deg' + ' '*1
    line += 'dGaussPosangle_deg' + ' '*1
    line += 'RMS_mJybm' + ' '*1
    line += 'Sky_mJybm' + ' '*1
    line += 'Sigma' + ' '*2
    line += 'mType' + ' '*1
    line += 'fArtefact' + ' '*1
    line += 'fCluster' + ' '*1
    line += 'fEdge' + ' '*1
    line += 'fHiNoise' + ' '*1
    line += 'fHi5Sig' + ' '*1
    line += 'fNearBright' + ' '*1
    line += 'fSmoothWeighting' + ' '*1
    line += 'fOverlap7sig' + ' '*1
    line += 'fOverlap5sig' + ' '*1
    if showSep:
        line += 'Sep_asec '
        line += 'Index Extra_text...'

    print line

    # Print the catalogue entries
    for e in rowLst:   
        line = ''
        line += '%s ' % e['sourceName']
        line += '%s ' % num2str(e['l_deg'],3,5,0,dlz,1)
        line += '%s ' % num2str(e['b_deg'],2,5,1,dlz,1)
        line += '%s ' % num2str(e['RAwcent_hr']*15.0, 3,5,0,dlz,1) 
        line += '%s ' % num2str(e['DecWcent_deg'],2,5,1,dlz,1)
        line += '%s ' % num2str(e['dRAabs_asec'],2,2,0,dlz,1) + ' '*3
        line += '%s ' % num2str(e['dDecAbs_asec'],2,2,0,dlz,1) + ' '*4
        line += '%s ' % num2str(e['mergePeakIntAbs_mJybm'],5,2,0,dlz,1) + ' '*2
        line += '%s ' % num2str(e['dMergePeakIntAbs_mJybm'],5,2,0,dlz,1) + ' '*3
        line += '%s ' % num2str(e['mergeFluxAbs_mJy'],5,2,0,dlz,1) 
        line += '%s ' % num2str(e['dMergeFluxAbs_mJy'],5,2,0,dlz,1) + ' '*1
        line += '%s ' % num2str(e['mergeAngscaleAbs_asec'],2,3,0,dlz,1) + ' '*7
        line += '%s ' % num2str(e['dMergeAngscaleAbs_asec'],2,3,0,dlz,1) + ' '*8
        if e['mergeAngscaleAbs_asec'] >= unresolved_asec:
            angDecon_asec = m.sqrt(e['mergeAngscaleAbs_asec']**2.0 - 2.25)
            line += '%s ' % num2str(angDecon_asec,2,1,0,dlz,1) + ' '*14
        else:
            line += ' '*19
        if e['mType']=='G':
            line += '%s ' % num2str(e['gaussMajorAxCor_asec'],2,3,0,dlz,1) + ' '*9
            line += '%s ' % num2str(e['dGaussMajorAxAbs_asec'],2,3,0,dlz,1) + ' '*10
            line += '%s ' % num2str(e['gaussMinorAxCor_asec'],2,3,0,dlz,1) + ' '*9
            line += '%s ' % num2str(e['dGaussMinorAxAbs_asec'],2,3,0,dlz,1) + ' '*10
            line += '%s ' % num2str(e['gaussPosangle_deg'],3,3,0,dlz,1) + ' '*10
            line += '%s ' % num2str(e['dGaussPosangle_deg'],2,3,0,dlz,1) + ' '*12
        else:
            line += ' '*103
        line += '%s ' % num2str(e['photMADFMsky_mJybm'],2,2,0,dlz,1) + ' '*4
        line += '%s ' % num2str(e['photMedianSky_mJybm'],1,3,1,dlz,1) + ' '*3
        line += '%s ' % num2str(e['photSigma'],4,1,0,dlz,1)
        line += '%s ' % e['mType'] + ' '*4
        if e['artifactProb'] == 'Unlikely':
            line += '0 ' + ' '*8
        elif e['artifactProb'] == 'Possibly':
            line += '1 ' + ' '*8
        else:
            line += '2 ' + ' '*8
        if e['fCluster']==7:
            line += '1 ' + ' '*7
        else:
            line += '0 ' + ' '*7
        if e['fEdge']==1:
            line += '1 ' + ' '*4
        else:
            line += '0 ' + ' '*4
        if e['fHiNoise']==1:
            line += '1 ' + ' '*7
        else:
            line += '0 ' + ' '*7
        if e['fHi5Sig']==1:
            line += '1 ' + ' '*6
        else:
            line += '0 ' + ' '*6
        if e['fNearBright']==1:
            line += '1 ' + ' '*10
        else:
            line += '0 ' + ' '*10
        if e['fSmoothWeighting']==1:
            line += '1 ' + ' '*15
        else:
            line += '0 ' + ' '*15
        if e['fOverlap7sig']>0:
            line += '1 ' + ' '*11
        else:
            line += '0 ' + ' '*11
        if e['fOverlap5sig']>0:
            line += '1 ' + ' '*11
        else:
            line += '0 ' + ' '*11

        if showSep:
            line += '%s ' % num2str(e['dist'] * 3600.0,6,1,0,dlz,1)
        if e.has_key('inCatDict'):
            line += '%s ' % num2str(e['inCatDict']['indx'],4,0,0,dlz,1)
            line += ' '.join(e['inCatDict']['extra'])

        print line


#-----------------------------------------------------------------------------#
# Serve the catalogue query results in CSV format                             #
#-----------------------------------------------------------------------------#
def print_cat_tab_csv(rowLst):
                                
    # Print the column headers
    #showCal = True;
    showSep, showCal, showCom, showBers, showWP = check_type()
    raExt = 'ra_new'
    decExt = 'dec_new'
    
    line = ''
    line += "Index,"
    line += "ra_deg,"
    line += "dec_deg,"
    line += "ra_hms,"
    line += "dec_deg,"
    line += "NVSS_St,"
    line += "NVSS_St_e,"
    line += "SUMSS_St,"
    line += "SUMSS_St_e,"
    line += "AT20G_St_20G,"
    line += "AT20G_St_20G_e"
    if showSep:
        line += ',Ang_Sep'
    if showCal:
        line += ',Cali_Name'
        raExt = 'ra_new_Cal'
        decExt = 'dec_new_Cal'
    if showCom:
        line += ',Alpha'
        line += ',Red_Chi_squ'
        line += ',SED_Fit_Name'
        raExt = 'ra_new_Com'
        decExt = 'dec_new_Com'
    if showBers:
        line += ',Kuhr_Name'
        line += ',ThreeC_Name'
        line += ',Class'
        line += ',Alpha'
        line += ',z'
        raExt = 'ra_new_B'
        decExt = 'dec_new_B'
    if showWP:
        line += ',WP85_Name'
        line += ',Alt_Name'
        line += ',WP85_Flux'
        line += ',Class'
        line += ',Alpha'
        line += ',z'
        raExt = 'ra_new_WP'
        decExt = 'dec_new_WP'
    print line

    # Print the catalogue entries
    i = 0;
    for e in rowLst:   
        i +=1;
        line = ''
        line += "%d," % i
        line += "%s," % (e[raExt])
        line += "%s," % (e[decExt])
        line += "%s," % deg2hms(e[raExt])
        line += "%s," % deg2dms(e[decExt], doSign=True)
        line += "%s," % e['N_St']
        line += "%s," % e['N_e_St']
        line += "%s," % e['tot_flux_dens']
        line += "%s," % e['tot_flux_dens_uncer']
        line += "%s," % e['S_20Ghz']
        line += "%s" % e['S_20Ghz_e']

        if showSep:
            line += ',%.2f' % (m.sqrt(float(e['dist'])) * 3600.0)
        if showCal:
            line += ',%s' % e['C_Name']
        if showCom:
            line += ',%.2f' % e['alpha']
            line += ',%.2f' % e['redchisq']
            line += ',%s' % e['NVSS_Name'] 
        if showBers:
            redshift= str(e['z']);
            if redshift == '0.0':
                redshift = ''
            line += ',%s' % e['Source_Name_1Jy']
            line += ',%s' % e['ThreeC_Name']
            line += ',%s' % e['Class']
            line += ',%.2f' % e['alpha']
            line += ',%s' % redshift
        if showWP:
            redshift = str(e['z']);
            if redshift == '0.0':
                redshift = ''
            line += ',%s' % e['IAUname']
            line += ',%s' % e['AltName']
            line += ',%s' % e['S11cm']
            line += ',%s' % e['Otype']
            line += ',%.2f' % e['Sp_Index']
            line += ',%s' % redshift

        print line
        
        
#-----------------------------------------------------------------------------#
# Serve the catalogue query results as KVIS ann file                          #
#-----------------------------------------------------------------------------#
def print_cat_tab_kvis(rowLst):

    print "COORD W"
    print "COLOUR green"
    for e in rowLst:
        print "ELLIPSE W",
        print '%f' % (e['RAwcent_hr']*15.0),
        print '%f' % e['DecWcent_deg'],
        print '%f' % (e['mergeAngscaleAbs_asec']/3600.0),
        print '%f' % (e['mergeAngscaleAbs_asec']/3600.0),
        print "0.0"


#-----------------------------------------------------------------------------#
# Serve the catalogue query results as AIPS star file                         #
#-----------------------------------------------------------------------------#
def print_cat_tab_aips(rowLst):
    
    for e in rowLst:
        print '%s' % deg2dms(e['RAwcent_hr'], delim=' '),
        print '%s' % deg2dms(e['DecWcent_deg'], delim=' '),
        print '%f' % (e['mergeAngscaleAbs_asec']),
        print '%f' % (e['mergeAngscaleAbs_asec']),
        print "0 3"

        
#-----------------------------------------------------------------------------#
# Serve the catalogue query results as DS9 region file                        #
#-----------------------------------------------------------------------------#
def print_cat_tab_ds9(rowLst):

    for e in rowLst:
        print "fk5;circle ",
        print '%s' % deg2dms(e['RAwcent_hr'], delim=':'),
        print '%s' % deg2dms(e['DecWcent_deg'], delim=':'),
        print '%f' % (e['mergeAngscaleAbs_asec'])


#-----------------------------------------------------------------------------#
# Parse the coordinate string input by the user                               #
#-----------------------------------------------------------------------------#
def parse_coord_string(coord_string,coord_system='Equatorial J2000'):
    ra_hr = None
    ra_deg = None
    dec_deg = None
    l_deg = None
    b_deg = None
    err = 0

    # Useful regular expressions
    re_decimal = re.compile('^(\d+(\.{1}\d*){0,1})')
    re_signed_decimal = re.compile('^([+|-]{0,1}\d+(\.{1}\d*){0,1})')
        
    # Collapse multiple spaces & multiple delimiters.
    #coord_string = coord_string.lstrip()               # leading spaces
    #coord_string = coord_string.rstrip()               # trailing spaces
    ## coord_string =  re.sub('\+',' +',coord_string)     # catch coords wo spaces
    ## coord_string =  re.sub('\-',' -',coord_string)     #
    ## coord_string = re.sub(',+',',',coord_string)       # multiple ','
    ## coord_string = re.sub('\s+',' ',coord_string)      # multiple spaces      
    ## coord_string = re.sub('\s*,\s*',',',coord_string)  # spaces+','+spaces = ,
    ## coord_string = re.sub('\s*:\s*',':',coord_string)  # spaces+':'+spaces = :
    ## coord_string = re.sub('\s*h\s*','h',coord_string)  # spaces+'h'+spaces = h
    ## coord_string = re.sub('\s*d\s*','d',coord_string)  # spaces+'d'+spaces = d
    ## coord_string = re.sub('\s*m\s*','m',coord_string)  # spaces+'m'+spaces = m
    ## coord_string = re.sub('\s*s\s*','s',coord_string)  # spaces+'s'+spaces = s
    ## coord_string = re.sub('h,','h',coord_string)       # 'h'+',' = h
    ## coord_string = re.sub('d,','d',coord_string)       # 'd'+',' = d
    ## coord_string = re.sub('m,','m',coord_string)       # 'm'+',' = m
    ## coord_string = re.sub('s,','s',coord_string)       # 's'+',' = s
    ## coord_string = re.sub('s$','',coord_string)        # 's'+',' = s
    
    # Split on allowed delimeters
    delim = re.compile('[,| |:|h|d|m|s]')
    coords = delim.split(coord_string)
    
    # Parse equatorial coordinates
    if (coord_system == 'Equatorial J2000'):
        
	# Check for decimal format
        if len(coords)==2:
            mch = re_decimal.match(coords[0])
            if mch:
                ra_deg = float(mch.group(1))
            else:
                ra_deg = None
            mch = re_signed_decimal.match(coords[1])
            if mch:
                dec_deg = float(mch.group(1))
            else:
                dec_deg = None

	# Or convert to decimal format
        elif len(coords)==6:
	    ra_hr   = dms2deg(','.join(coords[0:3]))
            ra_deg  = ra_hr * 15.0
	    dec_deg = dms2deg(','.join(coords[3:6]))

        # An error formatting
        else:
	    return [ra_deg, dec_deg, l_deg, b_deg, 1]
            	
	# Convert Equtorial to Galactic coordinates
        try:
            l_rad, b_rad = slalib.sla_eqgal(m.radians(ra_deg),
                                            m.radians(dec_deg))
            l_deg = m.degrees(l_rad)
            b_deg = m.degrees(b_rad)
        except:
            return [ra_deg, dec_deg, l_deg, b_deg, 1]

    # Parse Galactic coordinates
    elif coord_system == 'Galactic':

        # Check for only 2 arguments
        if len(coords)!=2: err=1

	# Check for correct decimal format
        mch = re_signed_decimal.match(coords[0])
        if mch:
            l_deg = float(mch.group(1))
        else:
            l_deg = None            
        mch = re_signed_decimal.match(coords[1])
        if mch:
            b_deg = float(mch.group(1))
        else:
            l_deg = None
	# Enforce the correct sign conventions
        if l_deg and b_deg:
            if (l_deg>=-360.0 and l_deg<0.0):
                l_deg = 360.0 + l_deg
            elif l_deg==360.0: 
                l_deg = 0.0
            elif (l_deg<-360.0 or l_deg>360.0):
                l_deg=None
            if (b_deg<-90.0 or b_deg>90.0):
                b_deg=None

	# Convert Galactic coordinates to Equ J2000 (FK5) coordinates
        try:
            ra_rad, dec_rad = slalib.sla_galeq(m.radians(l_deg),
                                               m.radians(b_deg))
            ra_deg = m.degrees(ra_rad)
            dec_deg = m.degrees(dec_rad)
        except:
            return [ra_deg, dec_deg, l_deg, b_deg, 1]

    # Coord system not recognised
    else:
        return [ra_deg, dec_deg, l_deg, b_deg, 1]
    
    # Final error check
    if (ra_deg == None or dec_deg == None): err=1
    if (l_deg == None or b_deg == None): err=1
    
    return [ra_deg, dec_deg, l_deg, b_deg, err]


#-----------------------------------------------------------------------------#
# Parse the file for coordinates                                              #
#-----------------------------------------------------------------------------#
def parse_coord_file_cone(fileHandle, coordSystem, maxLines=None):
    
    # Compile a few useful regular expressions
    comment = re.compile('#.*')
    spaces = re.compile('\s+')
    comma_and_spaces = re.compile(',\s+')

    # Process the file line-by-line
    validCounter = 0
    lineCounter = 0
    rowLst = []
    for line in fileHandle:
        entry = {}
        line = line.rstrip("\n\r")
        err = None
        
        # If not a comment 
        if not comment.match(line):

            # Clean up the syntax
            line = comment.sub('',line)
            line = line.strip()              # kill external whitespace
            line = spaces.sub(' ', line)     # shrink internal whitespace
            line = line.replace("'", '')     # kill quotes
            
            # Split entries into an array
            line = line.split()

            # Check for correct number of entries.
            if len(line)==0:
                continue
            lineCounter += 1
            entry['indx'] = lineCounter
            if not len(line)>1:
                err="less than 2 entries in line."

            # Parse the coordinates
            if not err:
                coordString = ', '.join(line[:2])
                [ra_deg, dec_deg, l_deg, b_deg, err] = \
                        parse_coord_string(coordString, coordSystem)
                if err:
                    err = 'coordinate parse error.'
                else:
                    entry['RA_hr'] = ra_deg / 15.0
                    entry['Dec_deg'] = dec_deg
                    entry['l_deg'] = l_deg
                    entry['b_deg'] = b_deg
                    entry['err'] = 0
                
            # Append further information
            if not err:
                entry['extra'] = line[2:]
            else:
                entry['err'] = err

            # Add the entry
            rowLst.append(entry)
            validCounter += 1
            if maxLines:
                if validCounter==maxLines:
                    break

    return rowLst


#-----------------------------------------------------------------------------#
# Parse the coordinate file used by the batch cutter                          #
#-----------------------------------------------------------------------------#
def parse_coord_file_cutter(fileName, coordSystem, size_amin, sizeMin_deg,
                            sizeMax_deg, outsystem, maxLines=None):
    
    coordParsedStr = ''
    inFileExt = '.fits'
    outFileExt = '_CORNISH_5GHz'

    # Parse the size variable
    size_deg, err = parse_radius(size_amin, 'arcmin', sizeMax_deg, sizeMin_deg)
    if err:
        return ''

    # Compile a few useful regular expressions
    comment = re.compile('#.*')
    spaces = re.compile('\s+')
    comma_and_spaces = re.compile(',\s+')

    # Process the file line-by-line
    f = open(fileName, 'r')
    validCounter = 0
    lineCounter = 0
    nameCheckDict = {}
    for line in f:
        lineCounter += 1
        line = line.rstrip("\n\r")
        err = None
        # If not a comment 
        if not comment.match(line):

            # Clean up the syntax
            line = comment.sub('',line)
            line = line.strip()              # kill external whitespace
            line = spaces.sub(' ', line)     # shrink internal whitespace
            line = line.replace("'", '')     # kill quotes
            
            # Split entries into an array
            line = line.split()

            # Check for correct number of entries.
            if len(line)==0:
                continue
            if not len(line)>1:
                err="less than 2 entries in line."

            # Parse the coordinates
            if not err:
                coordString = ', '.join(line[:2])
                ra_deg, dec_deg, l_deg, b_deg, err = \
                      parse_coord_string(coordString, coordSystem)
                if err:
                    err="failed to parse coordinates. Check 'input "+ \
                         "coordinate system' setting"
                
            # Parse the source name
            if not err:                    
                if len(line)>2:
                    name = re.sub('[\\\|*|?|&|!|@|\||/]', '', line[2])
                if len(line)==2 or name=='':
                    if b_deg>=0:
                        name = "G%03.3f+%03.3f" % (l_deg,b_deg,)
                    else:
                        name = "G" + "%04.3f" % l_deg + "%03.3f" % b_deg
                    name = str(lineCounter) + '_' + name
                
            # Avoid naming collisions
            if not err:
                if not nameCheckDict.has_key(name):
                    nameCheckDict[name] = 0
                else:
                    nameCheckDict[name]+=1
                if nameCheckDict[name]:
                    uniqueName = name+'_'+str(nameCheckDict[name])
                else:
                    uniqueName = name

            if outsystem == 'EQUATORIAL':
                x_deg = ra_deg
                y_deg = dec_deg
            else:
                x_deg = l_deg
                y_deg = b_deg

            # Record the parsed string or an error comment
            if err:
                coordParsedStr += "%d Error %s\n" % (lineCounter, err)
            else:
                coordParsedStr += "%d %s %s %s %.6f %.6f %.6f %s\n" \
                                  % (lineCounter, uniqueName, inFileExt,
                                     outFileExt, x_deg, y_deg, size_deg,
                                     outsystem)
                validCounter += 1
                if maxLines:
                    if validCounter==maxLines:
                        break
                    
    return coordParsedStr


#-----------------------------------------------------------------------------#
# Parse the radius value according to unit                                    #
#-----------------------------------------------------------------------------#
def parse_radius(size, unit, max_size_deg, min_size_deg=0.0): 
    err = 0
    size_deg = None
    size = size.lstrip()               # leading spaces
    size = size.rstrip()               # trailing spaces
    try:
        if unit == 'arcsec':
            size_deg = float(size)/3600.0            
        elif unit == 'arcmin':
            size_deg = float(size)/60.0
        else:
            size_deg = float(size)
    except:
        return [0.0, 1]
        
    # Size limits
    if size_deg > max_size_deg:
        size_deg = max_size_deg

    if size_deg < min_size_deg:
        size_deg = min_size_deg
    
    # Maximum size
    return size_deg, err


#-----------------------------------------------------------------------------#
# Map pixel coordinate to world coordinate on a single PNG image              #
#-----------------------------------------------------------------------------#
def get_clicked_coord(pngFile, fitsFile, xPixRaw, yPixRaw, offPixX=127,
                      offPixY=64, nPixXPNG=696, nPixYPNG=696):

    err = 0
    ra_deg = None
    dec_deg = None

    # Determine pixel clicked on the data
    png = im.open(pngFile)
    header = pf.getheader(fitsFile)
    scaleX = float(header['NAXIS1'])/float(nPixXPNG)
    scaleY = float(header['NAXIS1'])/float(nPixYPNG)
    sizeX = png.size[0]+1
    sizeY = png.size[1]+1
    xPixData = (xPixRaw - offPixX +1)*scaleX
    yPixData = (sizeY - yPixRaw - offPixY -1)*scaleY
    del png
    if (xPixData > header['NAXIS1'] or  yPixData > header['NAXIS2'] or
        xPixData <1 or  yPixData < 1):
        err = 1

    # Map to the world coordinate
    if not err:
        try:
            wcs = pw.WCS(header)
            coordsWorld = wcs.wcs_pix2sky([(xPixData,yPixData,0.0,0.0)], 1)
            ra_deg = coordsWorld[0][0]
            dec_deg = coordsWorld[0][1]
            del wcs, header
        except Exception:
            err = 2
            
    return (ra_deg,dec_deg,err)

#-----------------------------------------------------------------------------#
# Check what sort of query was requested when displaying in CSV               #
#-----------------------------------------------------------------------------#
def check_type():
    form = cgi.FieldStorage()
    showSep = False
    if (form.has_key("rad")):
        showSep = True;
    else:
        showSep = False;

    if form.has_key(""):
        if  form['errs'].value == '1':
            showErrs = True
    showCal = False;
    if form.has_key("Cal"):
        if form['Cal'].value == 'True':
            showCal = True
    if form['title'].value == 'Calibrator':
        showCal = True;
    
    showBers = False;
    if form.has_key("Bers"):
        if form['Bers'].value == 'True':
            showBers = True
    if form['title'].value == 'BERSCat':
        showBers = True;
        
    showCom = False;
    if form.has_key("compact"):
        if form['compact'].value == 'True':
            showCom = True
    if form['title'].value == 'Compact':
        showCom = True;
    
    showWP = False;
    if form.has_key("WP85"):
        if form['WP85'].value == 'True':
            showWP = True
    if form['title'].value == 'WP85':
        showWP = True;

    return (showSep, showCal, showCom, showBers, showWP)
