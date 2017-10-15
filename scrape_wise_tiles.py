#!/import/www/askap_targets/cgi-bin/python/bin/python
#=============================================================================#
# NAME:     scrape_wise_tiles.py                                              #
#                                                                             #
# PURPOSE:  Use the mechanize module to fetch the WISE tiles from the IPAC    #
#           server.                                                           #
#                                                                             #
# UPDATED:  mid-2014 by M. Glowacki (originally by C. Purcell)                #
#=============================================================================#
#/usr/bin/env python
# Coordinate system of returned images
returnSystem = 'FK5 - Equatorial J2000'
#returnSystem = 'GAL - Galactic' # Very very very slow!!!!!

# Import python modules
import os, sys, re, time
sys.path.insert(0, '/import/www/askap_targets/cgi-bin/python/lib/python2.7/site-packages/')
import numpy as np
import math as m
import mechanize
import cookielib
import urllib
from bs4 import BeautifulSoup
import html5lib
import cgitb
#import astropy as ast
cgitb.enable()

# Do not buffer output
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

# Needed for correct operation at USyd.
#os.environ['http_proxy'] = 'http://www-cache.usyd.edu.au:8080'

#status = 'no'
mechanize._sockettimeout._GLOBAL_DEFAULT_TIMEOUT = 360.0;

#-----------------------------------------------------------------------------#
# Main rountine                                                               #
#-----------------------------------------------------------------------------#
def getWise(ra, dec):
    #new - takes ra and dec input from the image_display files! 
    # Code to spoof a browser from:
    # http://stockrt.github.com/p/emulating-a-browser-in-python-with-mechanize/
    
    # Browser
    br = mechanize.Browser()
    status = 'no'
    # Cookie Jar
    cj = cookielib.LWPCookieJar()
    br.set_cookiejar(cj)
    
    # Browser options
    br.set_handle_equiv(True)
    #br.set_handle_gzip(True)
    br.set_handle_redirect(True)
    br.set_handle_referer(True)
    br.set_handle_robots(False)
    
    # Follows refresh 0 but not hangs on refresh > 0
    br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=5)
    
    # Want debugging messages?
    #br.set_debug_http(True)
    #br.set_debug_redirects(True)
    #br.set_debug_responses(True)

    # User-Agent (this is cheating, ok?)
    br.addheaders = [('User-agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]
    

    #-------------------------------------------------------------------------#

    # Login to the IPAC service
    #print "> Logging in to IPAC"
    r = br.open('http://hachi.ipac.caltech.edu:8080/montage/mosaicLogin.html')
    br.select_form(nr=0)
    br.form['user']=''
    br.form['passwd']=''
    #for obvious reasons not making the above public :V Won't work without an account. 
    br.submit()
    
    # Delete all complete jobs in the pipeline
    delete_all_jobs(br, True)

    l_deg = float(ra)
    d_deg = float(dec)
    coords = []
    coords.append((l_deg, d_deg))
    
    # Loop through the tiles centres
    posCount = 0
    completeTimeLst = []
    for lb in coords:

        #print "\n> Processing coordinates %s" % str(lb)
        jobIDlst  = []
        outFileDict = {}
        # Delete all complete jobs in the pipeline
        delete_all_jobs(br, True)
        
        # Submit four jobs in sequence - one for each WISE band
        bands = ['WISE 3.4 micron', 'WISE 4.6 micron',
                 'WISE 12 micron', 'WISE 22 micron']
        bandCode = ['W3.4', 'W4.6', 'W12', 'W22']
        
        for i in range(len(bands)):
            label = '%.3f_%.3f_%s' % (lb[0], lb[1], bandCode[i])
            #'%s' % bandCode[i]
            #'%.1f_%.1f_%s' % (lb[0], lb[1], bandCode[i])
            new = 'http://www.physics.usyd.edu.au/askap_targets/cgi-bin/'
            outFile = 'images/WISE_' + label + '.fits'
            if os.path.exists(outFile):
                #print "> Skipping existing file '%s' ..." % outFile
                status = 'exists'
            else:                
                jobID =  request_cutout(br, l_deg=lb[0], b_deg=lb[1],
                                        band=bands[i], label=label,
                                        returnSystem=returnSystem)
                outFileDict[jobID] = outFile
                jobIDlst.append(jobID)
        #print 'Job IDs: %s' % jobIDlst

        
        # Query the status of the jobs
        # Loop while some jobs are pending
        startTime = time.time()
        #print "> Checking job status ...",
        while True:
            allComplete = True
            try:
                 statusDict = check_job_status(br, jobIDlst)
            except Exception:
                #print "*"
                continue
            for jobID in jobIDlst:
                if not (statusDict[jobID]['status']=='COMPLETED'
                        or statusDict[jobID]['status']=='MISSING'
                        or statusDict[jobID]['status']=='ERROR'
                        or statusDict[jobID]['status']=='SERVERABORT'):
                    allComplete = False
                else:
                    try:
                        #trying to save the fits files more often - indivually, rather than waiting for one at a time
                        r = br.open(statusDict[jobID]['resultLink'])
                        if os.path.exists(outFile):
                            continue;
                        else:
                            for link in list(br.links(text_regex=re.compile("Mosaic"))):
                                fitsLnk = '/'.join(link.base_url.split('/')[:-1]) + \
                                    '/' + link.url
                                outFile = outFileDict[jobID]
#print "Saving the following URL to '%s'" % outFile
                    #print str(fitsLnk)
                                br.retrieve(fitsLnk, outFile)
                    except Exception:
                        continue;
                
            if allComplete:
                deltaTime = time.time() - startTime
                 #print '\nAll jobs complete in %.2f min with status:' % \
                #      (deltaTime/60.0)
                #for jobID in jobIDlst:
                #    print "ID: %s, Status: %s." % (jobID,
                #                                   statusDict[jobID]['status'])
                break
            else:
                #print '.',
                time.sleep(15)
                   
        # When all jobs are complete download the results of each band and
        # save to the disk
        nJobs = 0
        for jobID in jobIDlst:
            
            if statusDict[jobID]['status']!='COMPLETED':
                continue
            try:
                 
                r = br.open(statusDict[jobID]['resultLink'])
                for link in list(br.links(text_regex=re.compile("Mosaic"))):
                    fitsLnk = '/'.join(link.base_url.split('/')[:-1]) + \
                              '/' + link.url
                    outFile = outFileDict[jobID]
                    #print "Saving the following URL to '%s'" % outFile
                    #print str(fitsLnk)
                    if os.path.exists(outFile):
                        os.remove(outFile)
                    #image = urllib.URLopener()
                    br.retrieve(fitsLnk, outFile)
                    #print "... done."
                    nJobs += 1
            except Exception:
                #print "Failed to find result for job %s %s" % (jobID, fitsLnk)
                print '';
                
        
        # Clean up
        posCount +=1
        deltaTime = time.time() - startTime
        if nJobs>0:
            completeTimeLst.append(deltaTime)
            #print 'Position complete in %.2f min.' % (deltaTime/60.0)
            avgTime = np.median(completeTimeLst)
            remTime = (len(coords) - posCount) * avgTime / 3600.0
            #print 'Estimated time to list completion is %.2f hrs.' % remTime
        #delete_all_jobs(br, True)
        #return status

#-----------------------------------------------------------------------------#
# Delete all jobs running in the mosaic service                               #
#-----------------------------------------------------------------------------#
def delete_all_jobs(br, doAbort=False):
    
    #print "> Deleting existing jobs"
    abortLst = []
    deleteLst = []
    r = br.open('http://hachi.ipac.caltech.edu:8080/romeGetMessage?mainproglink=%2Fmontage')
    forms = [f for f in br.forms()]
    if forms==[]:
        #print 'No jobs in the queue'
        return
    
    controls = forms[0].controls
    for control in controls:
        if control.type=='checkbox':
            items = control.get_items()
            for item in items:
                if item.attrs['name']=='abortid':
                    abortLst.append(item.attrs['value'])
                if item.attrs['name']=='deleteid':
                    deleteLst.append(item.attrs['value'])
    if doAbort and abortLst!=[]:
        #print 'Aborting running jobs: ', abortLst
        forms[0]['abortid'] = abortLst
    if deleteLst!=[]:        
        #print 'Deleting complete jobs:', deleteLst
        forms[0]['deleteid'] = deleteLst
    br.select_form(nr=0)
    br.submit()
    

#-----------------------------------------------------------------------------#
# Submit a request for a cutout                                               #
#-----------------------------------------------------------------------------#
def request_cutout(br, l_deg, b_deg, band, size_deg=0.1,
                   returnSystem='GAL - Galactic',
                   label=''):

    #print "> Submitting job request ...",
    r = br.open('http://hachi.ipac.caltech.edu:8080/montage/index.html')
    forms = [f for f in br.forms()]
    forms[0]['locstr'] = '%f %f' % (l_deg, b_deg)
    forms[0]['band'] = [band]
    forms[0]['size'] = '%f' % size_deg
    forms[0]['resolution'] = ['1.375 arcsec (WISE Allsky Release)']
    #forms[0].set_value_by_label([returnSystem], name='resolution')
    forms[0]['cframe'] = ['FK5 - Equatorial J2000']
    forms[0]['label'] = label
    br.select_form(nr=0)
    br.submit()

    # Parse the response
    response = br.response().read().replace('\n', '')
    jobIDRe = re.compile('.*Request ID: (\d+).*')
    mch = jobIDRe.match(str(response))
    if mch:
        jobID = mch.group(1)
        #print "success"
    else:
        jobID = None
        #print "fail"
        
    return jobID
        

#-----------------------------------------------------------------------------#
# Check the status of a job in the server                                     #
#-----------------------------------------------------------------------------#
def check_job_status(br, jobIDlst=[]):

    statusDict = {}
    
    r = br.open('http://hachi.ipac.caltech.edu:8080/romeGetMessage?mainproglink=%2Fmontage')

    # Extract job information from 6th table
    html = br.response().read().replace('\n', '').replace('\r', '')
    soup = BeautifulSoup(html)
    tab = soup.findAll('table')
    if len(tab) < 6:
        #print "No jobs in queue"
        return statusDict

    # Loop through the rows and extract the information
    rows = tab[5].findAll('tr')
    
    for tr in rows[1:]:
        e = {}    # Row entry dict 
        # Separate the text entries and HTML for each cell
        colsTxt = tr.findAll(text=re.compile('[\S\w]'))
        cols = tr.findAll('td')
        ###breaks in below line###
        e['status'] = colsTxt[1].strip()
        e['label'] = colsTxt[2].strip()
        # Find the results link and extract the href attribute
        resultLink = None
        if len(cols[3].findAll('a'))==1:                    
            resultLink = cols[3].findAll('a')[0]['href']
            
        e['resultLink'] = resultLink
        # The return dictionary is indexed by the job number
        statusDict[colsTxt[0].strip()] = e
    # Filter the output for input job IDs
    
    if len(jobIDlst)>0:
        filtStatusDict = {}
        eBlank = {'status': u'MISSING', 'resultLink': None, 'label': None}
        for jobID in jobIDlst:
            if jobID in statusDict:
                filtStatusDict[jobID] = statusDict[jobID]
            else:
                filtStatusDict[jobID] = eBlank
        return filtStatusDict
    else:
        
        return statusDict
        
#=============================================================================#
#getWise()

