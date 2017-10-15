#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_html.py                                                      #
#                                                                             #
# PURPOSE:  Utility HTML formatting functions for serving CGI web pages.      #
#                                                                             #
# MODIFIED: 20-Mar-2013 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

# Import standard python modules
import os
import re
import sys
import shutil

import cgi
from cgi import escape
import Cookie
import cgitb
cgitb.enable()
import urllib


#-----------------------------------------------------------------------------#
# Format a number with leading spaces/zeros/signs.                            #
#-----------------------------------------------------------------------------#
def num2str(number,pre,post,dosign=False,doZeros=False,doPlus=False):
    lead=''
    formatCode=''
    if doZeros:
        lead = '%0'
    else:
        lead = '%'
    if dosign:
        if abs(number)/number>0.0:
            if doPlus:
                lead = '+' + lead
            formatCode = lead + str(pre+post+1) + ".%df" % post
        else:
            formatCode = lead + str(pre+post+2) + ".%df" % post
    else:
        formatCode = lead + str(pre+post+1) + ".%df" % post
    return "%s" % (formatCode) % number


#-----------------------------------------------------------------------------#
# Format a number with leading spaces/zeros/signs.                            #
#-----------------------------------------------------------------------------#
def str2cgi(string):
    return urllib.quote(string)


#-----------------------------------------------------------------------------#
# Print the HTML header for the page                                          #
#-----------------------------------------------------------------------------#
def print_page_head(stylesheet, title='', autoRefresh=0):
    print '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"',
    print '"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n'
    print '<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" '
    print 'xml:lang="en-US">\n'
    print '<head>\n'
    print '<title>%s</title>\n' % title
    print '<link rel="stylesheet" type="text/css" href="%s"/>\n' % stylesheet
    print '<meta http-equiv="Content-Type" content="text/html; ',
    print 'charset=iso-8859-1"/>'
    if autoRefresh:
        print '<META HTTP-EQUIV="REFRESH" CONTENT="%d">' % autoRefresh
    print '</head>\n'


#-----------------------------------------------------------------------------#
# Print a file in the current directory (i.e., include external html code)    #
#-----------------------------------------------------------------------------#
def include_file(filename):
    if os.path.exists(filename):
        FILEH = open(filename,'r')
        for line in FILEH:
            print line
        FILEH.close()


#-----------------------------------------------------------------------------#
# Convert HH/DD:MM:DD.D coordinates to decimal degrees.                       #
#-----------------------------------------------------------------------------#
def hms2deg(hms, store=':'):
    try:
        #delim=re.compile('[,| |:|h|d|m|s]')
        (h,m,s) = hms.split(store)
        h=float(h)
        m=float(m)
        s=float(s)
        #sign = 1
        #if d != 0.0:
        #    sign = int(d/(abs(d)))
        return (abs(h)*15+m/4.0+s/240.0)
    except Exception:
        return None

def dms2deg(dms, store=':'):
    try:
        #delim=re.compile('[,| |:|h|d|m|s]')
        (d,m,s) = dms.split(store)
        d=float(d)
        m=float(m)
        s=float(s)
        sign = 1
        if d != 0.0:
            sign = int(d/(abs(d)))
        return sign*(abs(d)+m/60.0+s/3600.0)
    except Exception:
        return None

#-----------------------------------------------------------------------------#
# Convert a float in degrees to 'hh mm ss' format                             #
#-----------------------------------------------------------------------------#
def deg2hms(deg, delim=':', nPlaces=2):

    try:
        angle = abs(deg)
        #sign=1
        #if angle!=0: sign = angle/deg
        
        # Calcuate the degrees, min and sec
        hh = int(angle/15)
        mm = int((angle-15*hh)*4)
        ss = (4*deg-60*hh-mm)*60

        # If rounding up to 60, carry to the next term
        if float("%05.2f" % ss) >=60.0:
           mm+=1.0
           ss = ss - 60.0
        if float("%02d" % mm) >=60.0:
           hh+=1.0
           mm = mm - 60.0
        if nPlaces> 0:
            formatCode = "%0" + "%s.%sf" % (str(2 + nPlaces + 1), str(nPlaces))
        else:
            formatCode = "%02.0f"
        formatCode = "%02d%s%02d%s" + formatCode
        return formatCode % (hh, delim, mm, delim, ss)
        
    except Exception:
        return None


#-----------------------------------------------------------------------------#
# Convert a float in degrees to 'dd mm ss' format                             #
#-----------------------------------------------------------------------------#
def deg2dms(deg, delim=':', doSign=False, nPlaces=2):

    try:
        angle = abs(deg)
        sign=1
        if angle!=0: sign = angle/deg
        
        # Calcuate the degrees, min and sec
        dd = int(angle)
        rmndr = 60.0*(angle - dd)
        mm = int(rmndr)
        ss = 60.0*(rmndr-mm)

        # If rounding up to 60, carry to the next term
        if float("%05.2f" % ss) >=60.0:
            mm+=1.0
            ss = ss - 60.0
        if float("%02d" % mm) >=60.0:
            dd+=1.0
            mm = mm -60.0
        if nPlaces> 0:
            formatCode = "%0" + "%s.%sf" % (str(2 + nPlaces + 1), str(nPlaces))
        else:
            formatCode = "%02.0f"
        if sign>0:
            if doSign:
                formatCode = "+%02d%s%02d%s" + formatCode
            else:
                formatCode = "%02d%s%02d%s" + formatCode
        else:
            if doSign:
                formatCode = "-%02d%s%02d%s" + formatCode
            else:
            	formatCode = "%02d%s%02d%s" + formatCode
        return formatCode % (dd, delim, mm, delim, ss)
        
    except Exception:
        return None

    
#-----------------------------------------------------------------------------#
# Store the results of a select statement into a list of dictionaries         #
#-----------------------------------------------------------------------------#
def select_into_namedlst(cursor, sql, args=[]):

    rowLst = []
    if args == []:
        cursor.execute(sql)
    else:
        cursor.execute(sql, tuple(args))
    columnNameLst = zip(*cursor.description)[0]
    rows = cursor.fetchall()
    for row in rows:
        e = {}
        for i in range(len(columnNameLst)):
            e[columnNameLst[i]] = row[i]
        rowLst.append(e)
    return rowLst


#-----------------------------------------------------------------------------#
# Convert a list of dictionaries into a dictionary of dictionaries            #
#-----------------------------------------------------------------------------#
def rowlst_to_rowdict(rowLst, keyName):
    rowDict = {}
    for e in rowLst:
        key = e.pop(keyName)
        rowDict[key] = e

    return rowDict

