#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     cat_query.py                                                      #
#                                                                             #
# PURPOSE:  Perform pre-set queries on the NVSS/SUMSS base  database.         #
#                                                                             #
# MODIFIED: some date by M. Glowacki (Originally by C. Purcell)               #
#                                                                             #
#=============================================================================#

#import
import os
import sys
import re
import math as m
#so my python set-up is used, NOT the uni's!
sys.path.append('/usr/lib/python2.7/dist-packages/')
sys.path.append('/usr/local/lib/python2.7/dist-packages/')
os.environ[ 'HOME' ] = '/tmp/'
import sqlite3
import cgi
from cgi import escape
import cgitb
cgitb.enable()

from util_html import *
from util_ASKAP import *
conn = sqlite3.connect('master.db')
cursor = conn.cursor()

def main():
    #global passing of whether the query is for a specific survey or not, e.g. showCal -> ATCA Calibrators
    global showSep, showCal, showBers, showWP, showCom, form

    #get the paths
    rootURL = open('env_rootURL.txt').read().strip()
    cgiURL = open('env_cgiURL.txt').read().strip()
    #opens the file, reads in the file,
    #strips in the case there was whitespace or something else
    
    #this function does all the cgi stripping for these variables
    strong, showSep, showErrs, showCal, showBers, showWP, showCom, queryTitle, retType = cgiStrip();
    
    try:
        sql, pageTitle = get_preset_query(queryTitle, strong)
    except:
        return
    #If there is an error, stop running after error page is served.
    #runs another function given the above within this file
    if sql is None:
        serve_error(rootURL, 'undefined_query', queryTitle)
        return
    #if there is nothing successfully returned to sql,#give a relevant server error.

    try:
	rowLst = select_into_namedlst(cursor, sql)
    except:
        serve_error(rootURL, 'bad_coords')
	cursor.close()
    	conn.close()
        return
    #list of rows using output of sql query, via that function.
    #function in util_html
    
    # Disconnect from the MySQL database
    cursor.close()
    conn.close()
    #not committing changes, naturally
    cc = 0;

    if len(rowLst)<1:
        serve_error(rootURL, 'zero_sources')
        return
    #if there's no results returned, display relevant message

    #otherwise, display results depending on the requested format (atm, HTML and CSV only)
    else:
        #if retType == 'ASCII':
        #    print "Content-type: text/plain\n"
        #    print_cat_tab_ascii(rowLst)
        if retType == 'CSV':
            print "Content-type: text/plain\n"
            print_cat_tab_csv(rowLst)
        #elif retType == 'KVIS':
        #    print "Content-type: text/plain\n"
        #    print_cat_tab_kvis(rowLst)
        #elif retType == 'DS9':
        #    print "Content-type: text/plain\n"
        #    print_cat_tab_ds9(rowLst)
        #elif retType == 'AIPS':
        #    print "Content-type: text/plain\n"
        #    print_cat_tab_aips(rowLst)
        else:
            serve_page_html(rootURL, cgiURL, rowLst, pageTitle, showErrs, showCal, showCom, showBers, showWP)
 #these are in util_ASKAP.py, i.e. util_CORNISH.py
    
#-----------------------------------------------------------------------------#
# Serve the catalogue query results as a web-page                             #
#-----------------------------------------------------------------------------#

def serve_page_html(rootURL, cgiURL, rowLst, pageTitle='Query Results', showErrs=False, showCal=False, showCom=False, showBers=False, showWP=False):
    #inputs are the root and cgi url, the rows obtained, and page title string
    # Print the page header
    print "Content-Type: text/html\n\n"
    print_page_head(rootURL + '/style.css', 'NVSS & SUMSS Catalogue Query')
    #style sheet largely same as Cormac's, may change in future
    ##include_file('includes/header_nonav.html')
    #if I want one, the header file will be different, currently do not have anyway

    # Print the title
    print '<h2><p style="text-align:center;">%s</h2><br>' % pageTitle
    print '''
    <div align="center">
    <h3>Search Results:</h3>
    <p>The query returned %d NVSS and SUMSS sources.</p>
    ''' % (len(rowLst))
    #changed above from CORNISH. len(rowLst) = number of results.
            
    # Print the HTML to serve the catalogue table - via function
    print_cat_html_alt(rootURL, cgiURL, rowLst, showSep, showErrs=showErrs, showCal=showCal, showCom=showCom, showBers=showBers, showWP=showWP)

    print '</div>'
    #Code for a hit tracker. Yay for javascript. 
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

    print '</body>\n'
    print '</html>\n'

#-----------------------------------------------------------------------------#
# Serve an error page with the relevant error message                         #
#-----------------------------------------------------------------------------#
def serve_error(rootURL, err, errStr=''):
    #This needs more work!
    # Print the page header
    print "Content-Type: text/html\n\n"
    print_page_head(rootURL + '/style.css', 'Catalogue Query')
    ##include_file('includes/header_nonav.html')
    #CHANGE!

    if err == 'blank_field':
        print '''
        <center><h3>No query title specified in the URL.</h3></center>
        <center><p>Please go back and try again.</p></center>
        '''
    elif err == 'undefined_query':
        print '''
        <center><h3>The requested catalogue query is not defined!</h3></center>
        <center><p>The query name in the URL was "title=%s". Please contact the
        system administrator and report this error, after ensuring the query title is correct.</p></center>
        ''' % errStr
    elif err == 'zero_sources':
        print '''
        <center><h3>No sources found in database for this query!</h3></center>
        '''
    elif err == 'bad_coords':
        print '''
        <center><h3>Please check your inputted coordinates or settings - wrong format or illegal characters detected!</h3></center>
        '''
    else:
        print """<center><h3>There's something nasty in the woodshed ...
        </h3></center>
       <center><p>An unknown error occured!</p></center>
       """

    print '</body>\n'
    print '</html>\n'

#-----------------------------------------------------------------------------#
# Choose the SQL and page title based on the query title                      #
#-----------------------------------------------------------------------------#
def get_preset_query(queryTitle, strong=False):
    global showCal, showBers, showWP, showCom
    valid = 0
    pageTitle = '';
    raExt = '';
    #raExt is null when just looking at the NVSS/SUMSS base table (columns of ra_new, dec_new)
    #the nme of the column changes because of difficulties crossmatching in sqlite - names must differ!
    #e.g. 'ra_new_Cal' for the RA in the table of ATCA Calibrator sources matched with NVSS/SUMSS table. 

    selSQL, fromSQL = get_base_cat_sql()
    #the base SQL query from this function
    #from util_CORNISH.py
    sql = selSQL + fromSQL;
    sql = selSQL + ', ra_new, dec_new ' + fromSQL
    sql = sql.replace('\n', ' ')
    #arranges the SQL query
    rootURL = open('env_rootURL.txt').read().strip()
    form = cgi.FieldStorage();    
    
    #Many 'if this is the query, this will be the SQL for said query' IF statements
    #Boring but it works. 

    #Query for returning objects within a range of RA and DEC positions
    if queryTitle == 'Pos_Range':
        valid = 1
        case = 2;
        #where case = the format of the given coordinates
        try:
            (rmin, rmax, dmin, dmax) = con_coords(case)
        #this function runs code for getting the ra, dec min/max values
        #and converts if given in hms/dms or otherwise

        #If users put the RA/DEC min/max values in the wrong order, fix that!
        #change later to a fancier/quicker min/max(form[].value) deal
            if float(rmin) > float(rmax):
                rmax = form['ra_min'].value
                rmin = form['ra_max'].value
                if float(dmin) > float(dmax):
                    dmax = form['dec_min'].value
                    dmin = form['dec_max'].value        

            posList = (float(rmin), float(rmax), float(dmin), float(dmax))
            pageTitle = "Radio Sources in region of sky [Ra: %0.2f to %0.2f degrees; Dec: %0.2f to %0.2f degrees]" % posList
            sql += ' WHERE ra_new > %f AND ra_new < %f AND dec_new > %f AND dec_new < %f' % posList
            sql = sql.replace('\n', ' ')
        except:
            serve_error(rootURL, 'bad_coords')
            return
    
    #Query for objects within a given radius of a specific position
    if queryTitle == 'Pos_Search':
        valid = 1
        case = 1;
        try:
            (ra, dec, rad) = con_coords(case);
            pageTitle = "Radio Sources within %0.2f arcminutes of (%0.2f, %0.2f)" % (float(rad), float(ra), float(dec))
            rad = str(float(rad)/60.0)
            sqlNew = query_catalogue_cone(ra, dec, rad)
            sql = sqlNew;
            sql = sql.replace('\n', ' ')
        #see function. Converts hms/dms to deg, etc.
        except:
            serve_error(rootURL, 'bad_coords')
            return
        #convert to degrees
        #function for thisthe SQL string is in a separate file
        

    #Query for the strongest flux objects
    if queryTitle == 'Strongest':
        try:
	    valid = 1;
            pageTitle = "Strongest Radio Sources";
	    form = cgi.FieldStorage();
	    order = form['order'].value
	    num = form['num'].value
	    sql = selSQL + ', ra_new, dec_new ' + fromSQL
	    #order will be NVSS or SUMSS - user decides to rank by NVSS or SUMSS
	    if order == 'NVSS':
	        sql += ' WHERE N_St BETWEEN 200 AND 1e8'
	        sql += ' ORDER BY N_St DESC LIMIT %f' % float(num);
	    if order == 'SUMSS':
	        sql += ' WHERE tot_flux_dens BETWEEN 200 AND 1e8'
	        sql += ' ORDER BY tot_flux_dens DESC LIMIT %f' % float(num);
	except:
            serve_error(rootURL, 'bad_coords')
            return
    
    #Query for the ATCA Calibrator Catalogue
    if queryTitle == 'Calibrator':
        valid = 1;
        raExt = '_Cal';
        try:
	    pageTitle = "Sources within ATCA Calibrator Catalogue"
            sql = selSQL + ', ra_new_Cal, dec_new_Cal, C_name ' + 'FROM Calibrator t1'#fromSQL
        #sql += ' INNER JOIN Calibrator'
        #sql += ' ON complete_all.ra_new = Calibrator.ra_com'
        #sql += ' AND complete_all.dec_new = Calibrator.dec_com'
        #changed because of decision to go to different ra/dec variables here
        # different table, else match elsewhere screwed over with sqlite!
            sql += ' WHERE ra_new_Cal BETWEEN 0 AND 360'
	except:
            serve_error(rootURL, 'bad_coords')
            return
    
    #Query for the sources in NVSS found by Joe Callingham to have |alpha| < 0.5 (i.e. likely compact sources)
    if queryTitle == 'Compact':
        valid = 1;
        raExt = '_Com';
	try:        
	    pageTitle = "Compact sources within NVSS (|alpha| < 0.5)"                
            sql = selSQL + ', NVSS_Name, redchisq, alpha, alpha_err, Flag, ra_new_Com, dec_new_Com ' + 'FROM Compact t1'#fromSQL
        #sql += ' INNER JOIN Compact'
        #sql += ' ON complete_all.ra_new = Compact.ra_new_com'
        #sql += ' AND complete_all.dec_new = Compact.dec_new_com'
            sql += ' WHERE ra_new_Com > 0'
        #this may seem silly, but this ensures that all _matched_ objects _only_ 
        #between the catalogue and NVSS/SUMSS objects are displayed
        #Later I may add in code for another ra column for objects matched or not matched, if this is desirable for the site
        #repeated in other queries

        #decide the flags:
            if form.has_key("suspect"):
                if form['suspect'].value == 'on':
                    sql += ' AND ra_new_Com > 0';
                else:
                    sql += ' AND Flag = "g"'
        
        #decide the reduced chisquare limit:
            if form.has_key("redchi"):
                try:
                    if float(form['redchi'].value) > 0:
                        sql += ' AND redchisq < %f' % float(form['redchi'].value)
                except Exception:
                    sql += ''
            sql = sql.replace('\n', ' ')
	except:
            serve_error(rootURL, 'bad_coords')
            return

    #Query for the Kuhr 1-Jansky Catalogue. Referred to as BERS (Bright Extragalactic Radio Sources)
    if queryTitle == 'BERSCat':
        valid = 1;
        raExt = '_B';
        try:
	    pageTitle = "Sources within the Bright Extragalactic Radio Sources (>1 Jy) Catalogue"
            sql = selSQL + ', Source_Name_1Jy, ThreeC_Name, Class, alpha, z, ra_new_B, dec_new_B ' + 'FROM Bers t1'#fromSQL;
        #sql += ' INNER JOIN Bers'
        #sql += ' ON complete_all.ra_new = Bers.ra_new_B'
        #sql += ' AND complete_all.dec_new = Bers.dec_new_B'
            sql += ' WHERE ra_new_B BETWEEN 0 AND 360'
	except:
            serve_error(rootURL, 'bad_coords')
            return

    #Query for the Wall & Peacock (1985) 2-Jansky catalogue
    if queryTitle == 'WP85':
        valid = 1;
        raExt = '_WP';
        try:
	    pageTitle = "Sources within the Wall and Peacock 2-Jy Radio Catalogue"
            sql = selSQL + ', IAUname, AltName, Otype, S11cm, Sp_Index, z, ra_new_WP, dec_new_WP' + ' FROM WP85 t1'
            sql += ' WHERE ra_new_WP BETWEEN 0 and 360'
	except:
            serve_error(rootURL, 'bad_coords')
            return
    #No longer used because the other way (a separate .py file for the Mixed query and display) was simplier.
    if queryTitle == 'Mixed':
        valid = 0;
	try:        
	    showCal = True;
            showBers = True
            pageTitle = "Sources within any additional table to NVSS/SUMSS/AT20G"
            defin = ' INNER'
        #defin = ' LEFT OUTER'
        #temp = selSQL.replace('N_St,', ' ')
            sql = selSQL + '''
, C_name, t1.ra_new_Cal, t1.dec_new_Cal, t2.Source_Name_1Jy, t2.ThreeC_Name, t2.Class, t2.alpha, t2.z 
''' + ' FROM Calibrator t1' 
            sql += ' %s JOIN Bers t2' % defin
            sql += ' ON t1.ra_new_Cal = t2.ra_new_B'
            sql += ' AND t1.dec_new_Cal = t2.dec_new_B'
	except:
            serve_error(rootURL, 'bad_coords')
            return
     #This concludes the queries - now for extra SQL for e.g. minimum flux limit

    sValue = 200.0 #in mJy; minimum value for strong sources
    #following is a way to impose a higher limit based on user input/selection
    
    if form.has_key("sValue"):
        try:
            sValue = float(form['sValue'].value)*1000.0
            #print sValue
        except Exception:
            sValue = 200.0;
        if sValue < 200.0:
            sValue = 200.0;

    #the strong query
    if queryTitle == 'strong':
        valid = 1
        #strong is a variable, turned on if this query is requested, or if the user wants strong sources
        strong = True
        try:
	    pageTitle = "Strong radio sources in NVSS and SUMSS"
            sql += ' WHERE ra_new > 0'
        #sql += ' WHERE (N_St BETWEEN 200 AND 1e8)'
        #sql += ' OR (tot_flux_dens BETWEEN 200 AND 1e8)'
        #sql += ' OR (S_20Ghz BETWEEN 200 AND 1e8)'
        #CONTINUED below within other loops, depending on various factors
	except:
            serve_error(rootURL, 'bad_coords')
            return

    #limits on declination:
    posLim = 'off'
    if form.has_key("poslim"):
        try:
            posLim = form['poslim'].value
        except Exception:
            posLim = 'off'
        try:
            decMax = float(form['decmax'].value)
        except Exception:
            decMax = 90;
        if decMax < -90:
            decMax = 90;

    if posLim == 'on':
        if (queryTitle == 'Strongest' and strong == True):
            sql += ''
        #exception is so the strong sql query is done properly; below statement in that case is given later
        else:
            try:
		sql += ' AND dec_new%s < %f' % (raExt,decMax)
	    except:
            	serve_error(rootURL, 'bad_coords')
            	return

    if strong == True:
        try:
	    if queryTitle != 'Strongest':
                sql += ' AND ((N_St BETWEEN %f AND 1e8)' % sValue
                sql += ' OR (tot_flux_dens BETWEEN %f AND 1e8)' % sValue
                sql += ' OR (S_20Ghz BETWEEN %f AND 1e8))' % sValue
            else:
                sql = selSQL + ', ra_new, dec_new ' + fromSQL
        #order will be NVSS or SUMSS - user decides to rank by NVSS or SUMSS
                if order == 'NVSS':
                    sql += ' WHERE N_St BETWEEN 200 AND 1e8'
                    sql += ' AND ((tot_flux_dens BETWEEN %f AND 1e8)' % sValue
                    sql += ' OR (S_20Ghz BETWEEN %f AND 1e8))' % sValue
                    if posLim == 'on':
                        sql += ' AND dec_new%s < %f' % (raExt,decMax)
                    sql += ' ORDER BY N_St DESC LIMIT %f' % float(num);
                if order == 'SUMSS':
                    sql += ' WHERE tot_flux_dens BETWEEN 200 AND 1e8'
                    sql += ' AND ((N_St BETWEEN %f AND 1e8)' % sValue
                    sql += ' OR (S_20Ghz BETWEEN %f AND 1e8))' % sValue
                    if posLim == 'on':
                        sql += ' AND dec_new%s < %f' % (raExt,decMax)
                    sql += ' ORDER BY tot_flux_dens DESC LIMIT %f' % float(num);
	except:
            serve_error(rootURL, 'bad_coords')
            return
    
    #Last bit - ORDER BY SQL statements depending on the query and requested sorting by the user
    #needs to be at the end here
    try:
	if queryTitle == 'Pos_Search':
            sql += ' ORDER BY dist';
        if queryTitle != 'Pos_Search' or queryTitle != 'Strongest':
            if form.has_key("sort"):
                sort = form['sort'].value
                sql += ' ORDER BY %s' % (sort)#, raExt)
                if sort == 'N_St' or sort == 'tot_flux_dens':
                    sql += ' '##DESC'
                #ordering in descending order doesn't work, lists N/A values first... troublesome
                #fix later? this still works. 
    except:
        serve_error(rootURL, 'bad_coords')
        return

    if valid==0:
        return None, pageTitle
    #if nothing worked, nothing returned, get error later too
    sql = sql.replace('\n', ' ')

    return sql, pageTitle
#otherwise, return the sql to execute and the relevant page title.

#=============================================================================#
def con_coords(case=1):
    #converts coords into a usable form [degrees] (the downsides of various hms formats...)
    #(ideally) works if the layout is e.g. 00:00:00.00, or 00h00m00.00, 00 00 00.00...
    #case 1 is for searching for objects within a radius
    if case==1:
        form = cgi.FieldStorage()
        if not (form.has_key("ra") and form.has_key("dec") \
                and form.has_key("rad") and form.has_key("title")):
            serve_error(rootURL, 'blank_field')
            return None, pageTitle
        
        ra = form['ra'].value
        dec = form['dec'].value
        ra = urllib.unquote(ra)
        dec = urllib.unquote(dec)
        rad = form['rad'].value
        try:
            if ':' in ra:
                store = ':'
                ra = hms2deg(ra,store)
            #converts the hms format to degrees, in util files.
            #made another one for dms2deg as well (used below)
            elif ' ' in ra:
                store = ' '
                ra = hms2deg(ra, store)
            elif 'm' in ra or 'M' in ra:
                store = ' '
                ra = re.sub('[sS]', '', ra)
            #in case there is a s or S at the end, so not to confuse the code            
                ra = re.sub('[a-zA-Z\n\,]', ' ', ra)
                ra = hms2deg(ra,store)
                if ':' in dec:
                    store = ':'
                    dec = dms2deg(dec, store)            
                elif ' ' in dec:
                    store = ' '
                    dec = dms2deg(dec, store)
                elif 'm' in dec or 'M' in dec:
                    store = ' '
                    dec = re.sub('[sS]', '', dec)
            #in case there is a s or S at the end, so not to confuse the code            
                    dec = re.sub('[a-zA-Z\n\,]', ' ', dec)
                    dec = dms2deg(dec,store)
            return (ra, dec, rad)
        except:
            serve_error(rootURL, 'bad_coords')
    #case 2 is for returning objects within a range of positions (so within a rectangle)
    if case == 2:
        form = cgi.FieldStorage()
        if not (form.has_key("ra_min") and form.has_key("ra_max") \
                and form.has_key("dec_min") and form.has_key("dec_max") \
                and form.has_key("title")):
            serve_error(rootURL, 'blank_field')
            return None, pageTitle
        
        rmin = form['ra_min'].value
        rmax = form['ra_max'].value
        dmin = form['dec_min'].value
        dmax = form['dec_max'].value
        rmin = urllib.unquote(rmin)
        dmin = urllib.unquote(dmin)
        rmax = urllib.unquote(rmax)
        dmax = urllib.unquote(dmax)
        if ':' in rmin:
            store = ':'
            rmin = hms2deg(rmin,store)
        elif ' ' in rmin:
            store = ' '
            rmin = hms2deg(rmin, store)
        elif 'm' in rmin or 'M' in rmin:
            store = ' '
            rmin = re.sub('[sS]', '', rmin)
            #in case there is a s or S at the end, so not to confuse the code            
            rmin = re.sub('[a-zA-Z\n\,]', ' ', rmin)
            rmin = hms2deg(rmin,store)
        if ':' in rmax:
            store = ':'
            rmax = hms2deg(rmax,store)
        elif ' ' in rmax:
            store = ' '
            rmax = hms2deg(rmax, store)
        elif 'm' in rmax or 'M' in rmax:
            store = ' '
            rmax = re.sub('[sS]', '', rmax)
            #in case there is a s or S at the end, so not to confuse the code            
            rmax = re.sub('[a-zA-Z\n\,]', ' ', rmax)
            rmax = hms2deg(rmax,store)
            
        if ':' in dmin:
            store = ':'
            dmin = dms2deg(dmin, store)            
        elif ' ' in dmin:
            store = ' '
            dmin = dms2deg(dmin, store)
        elif 'm' in dmin or 'M' in dmin:
            store = ' '
            dmin = re.sub('[sS]', '', dmin)
            #in case there is a s or S at the end, so not to confuse the code            
            dmin = re.sub('[a-zA-Z\n\,]', ' ', dmin)
            dmin = dms2deg(dmin,store)
        if ':' in dmax:
            store = ':'
            dmax = dms2deg(dmax, store)            
        elif ' ' in dmax:
            store = ' '
            dmax = dms2deg(dmax, store)
        elif 'm' in dmax or 'M' in dmax:
            store = ' '
            dmax = re.sub('[sS]', '', dmax)
            #in case there is a s or S at the end, so not to confuse the code            
            dmax = re.sub('[a-zA-Z\n\,]', ' ', dmax)
            dmax = dms2deg(dmax,store)
        return (rmin, rmax, dmin, dmax)

def cgiStrip():
    # Read the variables from the URL and assigns them to variables if they are present
    form = cgi.FieldStorage()
    if not form.has_key("title"):
        serve_error(rootURL, 'blank_field')
        return
    #if the form doesn't have this key, then return a server error
    #this calls a function that displays a relevant error message
    #depending on the second variable (string, e.g. 'blank_field').

    queryTitle = form['title'].value
    #assigns the form's input to the variable queryTitle

    strong = False;
    if form.has_key("strong"):
        if form['strong'].value == 'on':
            strong = True;
    #this structure repeats itself
    showErrs = False
    if form.has_key("errs"):
        if  form['errs'].value == '1':
            showErrs = True
    #same as above, just with showing errors!
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

    retType = "HTML";
    if form.has_key("type"):
        retType = form['type'].value

    #return them all
    return  strong, showSep, showErrs, showCal, showBers, showWP, showCom, queryTitle, retType


main()


    

    
