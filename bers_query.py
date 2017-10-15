#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     bers_query.py                                                #
#                                                                             #
# PURPOSE:  Perform pre-set queries on the catalogue database.        #
#                                                                             #
# MODIFIED: 2014 by M. Glowacki (Originally by C. Purcell)               #
#                                                                             #
#=============================================================================#
#/import/www/askap_targets/cgi-bin/python/bin/python
import os
import sys
import re
import math as m
sys.path.append('/usr/lib/python2.7/dist-packages/')
sys.path.append('/usr/local/lib/python2.7/dist-packages/')
sys.path.append('/import/www/askap_targets/cgi-bin/python/lib/python2.7/dist-packages/')
sys.path.append('/import/www/askap_targets/cgi-bin/python/lib/python2.7/dist-packages/numpy')
sys.path.append('/import/www/askap_targets/cgi-bin/python/lib/python2.7/site-packages/numpy')
#sys.path.insert(0, '/import/www/askap_targets/cgi-bin/python/lib/python2.7/site-packages/')
#sys.path.insert(0, '/import/www/askap_targets/cgi-bin/python/lib/python2.7/dist-packages/')
#sys.path.append(os.path.split(sys.argv[0])[0])
os.environ[ 'HOME' ] = '/tmp/'
import sqlite3
import cgi
from cgi import escape
import cgitb
cgitb.enable()
import urllib
#from astropy.coordinates import ICRS, Galactic
#from astropy import units as u

from util_html import *
from util_ASKAP import *

conn = sqlite3.connect('master.db')
cursor = conn.cursor()

#top level function

def main():
    global sourceName

    #get the paths
    rootURL = open('env_rootURL.txt').read().strip()
    cgiURL = open('env_cgiURL.txt').read().strip()
    #these are files with urls that may need changing, see to it later
    #opens the file, reads in the file,
    #strips in the case there was whitespace or something else

    # Read the variables from the URL
    form = cgi.FieldStorage()
    if not form.has_key("title"):
        serve_error(rootURL, 'blank_field')
        return
    #if the form doesn't have this key, then return a server error
    #this calls a function that displays a relevant error message
    #depending on the second variable (string, e.g. 'blank_field').

    queryTitle = form['title'].value
    #assigns the form's input to the variable queryTitle
    sourceName = form['source'].value
    if ' ' in sourceName:
        sourceName = sourceName[0:4]+'+'+sourceName[5:len(sourceName)];
    
    try:
	sql, pageTitle = get_preset_query(queryTitle)
    except:
	return
    #runs another function given the above within this file
    if sql is None:
        serve_error(rootURL, 'undefined_query', queryTitle)
        return
    #if there is nothing successfully returned to sql
    #give a relevant server error.

    try:
	rowLst = select_into_namedlst(cursor, sql)
	cursor.close()
    	conn.close()
    except:
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
    #otherwise, display results depending on the requested format, e.g. HTML
    else:
        serve_page_htmly(rootURL, cgiURL, rowLst, pageTitle)
    
#-----------------------------------------------------------------------------#
# Serve the catalogue query results as a web-page                             #
#-----------------------------------------------------------------------------#

def serve_page_htmly(rootURL, cgiURL, rowLst, pageTitle='Query Results', showErrs=False, showCal=False, showCom=False, showBers=False):

    #inputs are the root and cgi url, the rows obtained, and page title string
    # Print the page header
    print "Content-Type: text/html\n\n"
    print_page_head(rootURL + '/style.css', 'NVSS & SUMSS Catalogue Query')
    #CHANGE for own style sheet later - unimportant now
##    include_file('includes/header_nonav.html')
    #CHANGE/CHECK - header file will be different, currently do not have anyway

    # Print the title
    print '<h2><p style="text-align:center;">%s</h2><br>' % pageTitle
    print '''
    <div align="center">
    <h3>Search Results:</h3>
    <p>The query returned %d fluxes for this object from the BERS Catalogue.</p>
    ''' % (len(rowLst))
    #changed above from CORNISH. len(rowLst) = number of results.
    # Print the HTML to serve the catalogue table - all here for sake of ease
    
    print '''
    <p>The row containing the maximum flux has been bolded.
    <a href="%sbers_query.py?title=BERSFlux&source=%s&sort=S_bers">Click to sort by flux.</a></p>
    ''' % (cgiURL, sourceName)
    print '''
    <table class="fancy" style="margin-left:auto; margin-right:auto;
    text-align:center;">
    <tbody>
    <tr>
    <th>No.</th>
    <th>Source Name</th>
    <th>Frequency (MHz)</th>
    <th>Flux (Jy)</th>
    <th>Flux error (Jy)</th>
    '''
    # Print the source entries
    i = 0
    maxFlux = 0
    for f in rowLst:
        if float(f['S_bers']) > float(maxFlux):
            maxFlux = f['S_bers']
    for e in rowLst:
        ## if i == 0:
        ##     i += 1;
        ##     continue;
        i+=1
        print '</tr>'
        if float(e['S_bers']) == float(maxFlux):
            print '<td><b>%d</b></td>' % i
            try:
                print '<td><b>%s</b></td>' % e['Source_Name_F']
                print '<td><b>%d</b></td>' % int(e['Freq'])
                print '<td><b>%.2f</b></td>' % float(e['S_bers'])
                print '<td><b>%.2f</b></td>' % float(e['e_S_bers'])
            except Exception:
                print '<td></td>'
                print '<td></td>'
                print '<td></td>'
                print '<td></td>'
                    
        else:
            print '<td>%d</td>' % i
            try:
                print '<td>%s</td>' % e['Source_Name_F']
                print '<td>%d</td>' % int(e['Freq'])
                print '<td>%.2f</td>' % float(e['S_bers'])
                print '<td>%.2f</td>' % float(e['e_S_bers'])
            except Exception:
                print '<td></td>'
                print '<td></td>'
                print '<td></td>'
                print '<td></td>'
        print '</tr>'

    print '''
    </tbody>
    </table>
    </br></br>
    '''
    print '</div>'

    # Print the page footer - CHANGE/CHECK
    ## include_file('includes/footer.html')
    print '</body>\n'
    print '</html>\n'

#-----------------------------------------------------------------------------#
# Serve an error page with the relevant error message                         #
#-----------------------------------------------------------------------------#
def serve_error(rootURL, err, errStr=''):

    # Print the page header
    print "Content-Type: text/html\n\n"
    print_page_head(rootURL + '/style.css', 'CORNISH Catalogue Query')
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
        system administrator and report this error.</p></center>
        ''' % errStr
    elif err == 'zero_sources':
        print '''
        <center><h3>No matches found in database for source %s!</h3></center>
        ''' % sourceName
    elif err == 'bad_coords':
        print '''
        <center><h3>Please check your inputted coordinates or settings - wrong format or illegal characters detected!</h3></center>
        '''
    else:
        print """<center><h3>There's something nasty in the woodshed ...
        </h3></center>
       <center><p>An unknown error occured!</p></center>
       """

    # Print the page footer - CHANGE
   ## include_file('includes/footer.html')
    print '</body>\n'
    print '</html>\n'

#-----------------------------------------------------------------------------#
# Choose the SQL and page title based on the query title                      #
#-----------------------------------------------------------------------------#
def get_preset_query(queryTitle, strong=False):
    global sourceName

    valid = 0
    #boolean as a 0/1 case?
    pageTitle = ''
    #empty for now

    #arranges the SQL query
    rootURL = open('env_rootURL.txt').read().strip()
    form = cgi.FieldStorage();
    
    try:
	if queryTitle == 'BERSFlux':
	    valid = 1;
	    pageTitle = "Fluxes for source %s within the Bright Extragalactic Radio Sources (>1 Jy) Catalogue" % sourceName
            sql = 'SELECT Source_Name_F, Freq, S_bers, e_S_bers';
            sql += ' FROM Bers_vir'
        #sql += ' INNER JOIN Bers_vir'
        #sql += ' ON Bers.Source_Name_1Jy = Bers_vir.Source_Name_F'
            sql += ''' WHERE Source_Name_F MATCH '"%s"' ''' % sourceName

        if queryTitle != 'Pos_Search' or queryTitle != 'Strongest':
            if form.has_key("sort"):
            	sort = form['sort'].value
            	sql += ' ORDER BY %s DESC' % sort
            	if sort == 'N_St' or sort == 'tot_flux_dens':
                    sql += ' '##DESC'
                #ordering in descending order doesn't work, lists N/A values first... troublesome
                #fix later? this still works. 
    except:
        serve_error(rootURL, 'bad_coords')
        return
    #additional ordering for this query.
    if valid==0:
        return None, pageTitle
    #if nothing worked, nothing returned, get error later too
    sql = sql.replace('\n', ' ')

    return sql, pageTitle
#otherwise, return the sql to execute and the relevant page title.


main()


    

    
