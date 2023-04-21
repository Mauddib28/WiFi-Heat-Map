#!/usr/bin/env python2

"""
Copyright (c) 2016, Bliksem Labs B.V.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this 
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
   this list of conditions and the following disclaimer in the documentation 
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

###
# Imports
###

# Imports of libraries
import sys
import numpy
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy as np

# Attempt to import known versions of etree
try:
  from lxml import etree
except ImportError:
  try:
    # Python 2.5
    import xml.etree.cElementTree as etree
  except ImportError:
    try:
      # Python 2.5
      import xml.etree.ElementTree as etree
    except ImportError:
      try:
        # normal cElementTree install
        import cElementTree as etree
      except ImportError:
        try:
          # normal ElementTree install
          import elementtree.ElementTree as etree
        except ImportError:
          print("Failed to import ElementTree from any known place")

###
# Globals
###
debugBit = 0

###
# Functions
###

# Process the Kismet GPSXML into columns.

def parse_xml(filename):
	tree = etree.parse(open(filename, 'rb'))

	ts = []
	bssid = []
	signal = []
	lat = []
	lon = []
	walked_lon = []
	walked_lat = []

	for z in tree.findall('.//gps-point'):
		# A lon/lat filter might be applied here
		# if float(z.get('lon')) < 3.942:
		#	continue

		if z.get('bssid') == 'GP:SD:TR:AC:KL:OG':
			walked_lon.append(float(z.get('lon')))
			walked_lat.append(float(z.get('lat')))

		elif z.get('signal_dbm') is not None:
			bssid.append(z.get('bssid'))
			ts.append(int(z.get('time-sec')))
			lat.append(float(z.get('lat')))
			lon.append(float(z.get('lon')))
			signal.append(int(z.get('signal_dbm')))

	return (ts, bssid, signal, lat, lon, walked_lon, walked_lat,)


# Draw parsed data on a surface

def draw_data(ts, bssid, signal, lat, lon, walked_lon, walked_lat):

	# We create a grid of 1000x1000
	grid_x, grid_y = numpy.mgrid[min(walked_lon):max(walked_lon):1000j, min(walked_lat):max(walked_lat):1000j]

	# We want to draw all unique APs
	bssids = list(set(bssid))

	# For each BSSID...
	for s in bssids:
		points_lon = []
		points_lat = []
		values = []
		h = []
		
		# Apply all points on an intermediate surface
		# so we can distinct points where we were, without reception
		for i in range(0, len(bssid)):
			if bssid[i] == s:
				hc = hash((lon[i], lat[i]))
				if hc not in h:
					points_lon.append(lon[i])
					points_lat.append(lat[i])
					values.append(float(signal[i]))
					h.append(hash((lon[i], lat[i])))

		# Optional: apply -100dBm where we don't have gathered data
		for i in range(0, len(walked_lon)):
			hc = hash((walked_lon[i], walked_lat[i]))
			if hc not in h:
				points_lon.append(lon[i])
				points_lat.append(lat[i])
				values.append(float(-100))
				h.append(hash((walked_lon[i], walked_lat[i])))

		# Interpolate the data
		grid = griddata((points_lon, points_lat), numpy.array(values), (grid_x, grid_y), method='cubic')

		# Store the bitmap in the current folder.
		plt.show()
		plt.imsave('%s.png' % (s), grid.T)

		# Calculate the World File for use in Qgis
		a = ((max(walked_lon)-min(walked_lon))/1000)
		b = 0
		c = 0
		d = ((max(walked_lat)-min(walked_lat))/1000)
		e = min(walked_lon)
		f = min(walked_lat)

		# Write the World File
		open('%s.pngw' % (s), 'w').write('%.16f\n%d\n%d\n%.16f\n%.16f\n%.16f' % (a, b, c, d, e, f,))

# NOTE: The 'gpsxml' files do NOT contain the SSID names associated to the 'bssid' values.
#   -> Will need to read in the 'netxml' file to establish names of wireless APs
#   - Should be able to follow the same structure of the 'parse_xml()' function written earlier (by original code author)
# Borrowing code from: https://github.com/kurankat/kis2kml/blob/master/kis2kml.py

# Function for parsing the netxml file
#def parse_netxml( netxmlFile ):
    # NOTE: Searching for the following fields 'SSID' and 'BSSID' that exist within 'wireless-network' enteries
    #   -> Additionaly, there are separate enteries for 'wireless-client's (use for additional tracking?)
    #   - 'type' field under the 'wireless-network' tag says what kind of wireless tranmission was seen
    #       - Ex: 'probe' means client searching for AP(s), 'infrastructure' means an AP letting others know of its location
    #       - NOTE: Search to couple this with HackRF(??)

# Open xml file and load it into a etree.ElementTree object
# Create a list of tree nodes that contain infrastructure networks
# Parse xml into a list of dictionaries with all important network data
def load_nets_from_xml(xfile):
    global total_discovered, runtime
    netnodes = []
    netlist_dicts = []
    clientlist = []

    print("Reading network information from %s\n" %xfile)

    # Open Kismet .netxml file and load into list of nodes
    try:
        with open(xfile, 'rt') as kismet:
            try:
                tree = etree.parse(kismet)
            except etree.ParseError as xmlerr:
                print("\n*** ERROR ***  Problem parsing input file.")
                print("               Is it a Kismet netxml file?")
                print("               Python says: %s\n" % xmlerr)
                usage()
                sys.exit(2)
    except IOError as ioerr:
        print("\n*** ERROR ***  Cannot read input file. Does it exist?")
        print("\tPython says: %s\n" % ioerr)
        usage()
        sys.exit(2)

    for node in tree.iter('detection-run'):
        runtime = node.attrib.get('start-time')

    # For checking entries in a database.  Un-needed at this time
    #if runtime_exists():
    #    print("This detection run (%s) has already been imported" %runtime)
    #    sys.exit()

    netnodes = pop_xml_netlist(tree)
    # For each wireless network node, create a dictionary, and append it to
    # a list of network dictionaries
    for node in netnodes:
        netlist_dicts.append(populate_net_dict(node))
        populate_client_list(node, clientlist)
    total_discovered = len(netnodes)
    print("")

    return netlist_dicts, clientlist
    # netlist_dicts - Each top level item is a different wireless networks (wn_num)
    #   -> 'essid' - Name of wireless AP point
    #   -> 'bssid' - BSSID for the wireless AP point (NOTE: can be different for repeaters)

# Function to return a list of eTree nodes. Takes the whole XML tree as the
# argument and returns a list[] of nodes containing only infrastructire networks
def pop_xml_netlist(whole_tree):
    nodelist = []
    for node in whole_tree.findall('.//wireless-network'):
        if (node.attrib.get('type') == 'infrastructure'):
            nodelist.append(node)
    if len(nodelist) == 0:
        print ("\n+++ WARNING +++  "
               "There don't seem to be any wireless networks in your "
               "input file\n")
        usage()
        sys.exit()
    return nodelist

# Create a list of clients. Each client is a list of the router bssid,
# the client MAC and the client max_signal
def populate_client_list(wireless_node, client_list):

    for lev1 in wireless_node:
        if lev1.tag == 'BSSID':
            bssid = lev1.text
        if lev1.tag == 'wireless-client':
            cldata = []
            if lev1.attrib['type'] != 'fromds':
                cldata.append(bssid)
                for clientinfo in lev1:
                    if clientinfo.tag == 'client-mac':
                        cldata.append(clientinfo.text)
                    if clientinfo.tag == 'snr-info':
                        for snr in clientinfo:
                            if snr.tag == 'max_signal_dbm':
                                cldata.append(snr.text)
            if len(cldata) > 0:
                client_list.append(cldata)

# Populate values of network dictionary from xml node
def populate_net_dict(wireless_node):
    wn = make_net_dict()
    wn['wn_num'] = wireless_node.attrib['number']
    wn['first_seen'] = wireless_node.attrib['first-time']
    wn['last_seen'] = wireless_node.attrib['last-time']
    wn['placeholder_encryption'] = []
    wn['clients'] = []

    # Iterate through first-level nodes and fill values into empty dictionary
    for lev1 in wireless_node:

        # Append the MAC address of clients as a list
        if lev1.tag == 'wireless-client':
            for cinfo in lev1:
                for ctag in cinfo.iter('client-mac'):
                    wn['clients'].append(ctag.text)

        # Loop through second-level nodes in SSID
        if lev1.tag == 'SSID':
            for ssid_info in lev1:
                # assign multiple ecryption fields to a temporary list
                for e in ssid_info.iter('encryption'):
                    wn['placeholder_encryption'].append(e.text)

                if ssid_info.tag == 'type':
                    wn['ssid_type'] = ssid_info.text
                if ssid_info.tag == 'max-rate':
                    wn['max_speed'] = ssid_info.text
                if ssid_info.tag == 'packets':
                    wn['packets'] = ssid_info.text
                if ssid_info.tag == 'beaconrate':
                    wn['beaconrate'] = ssid_info.text
                if ssid_info.tag == 'wps':
                    wn['wps'] = ssid_info.text
                if ssid_info.tag == 'wps-manuf':
                    wn['wps_manuf'] = ssid_info.text
                if ssid_info.tag == 'dev-name':
                    wn['dev_name'] = ssid_info.text
                if ssid_info.tag == 'model-name':
                    wn['model_name'] = ssid_info.text
                if ssid_info.tag == 'model-num':
                    wn['model_num'] = ssid_info.text
                if ssid_info.tag == 'wpa-version':
                    wn['ssid_wpa_version'] = ssid_info.text
                if ssid_info.tag == 'essid':
                    if ssid_info.attrib['cloaked'] == 'true':
                        wn['essid'] = ""
                    else: # Replace some characters that cause problems in KML
                        tempessid = ssid_info.text
                        wn['essid'] = tempessid.replace('&', '').replace('<', \
                            '').replace('>', '')
                    wn['cloaked'] = ssid_info.attrib['cloaked']

        if lev1.tag == 'BSSID':
            wn['bssid'] = lev1.text
        if lev1.tag == 'manuf':
            wn['manuf'] = lev1.text
        if lev1.tag == 'channel':
            wn['channel'] = lev1.text
        if lev1.tag == 'maxseenrate':
            wn['maxseenrate'] = lev1.text

        # Loop through snr information
        if lev1.tag == 'snr-info':
            for snr_info in lev1:
                if snr_info.tag == 'max_signal_dbm':
                    wn['max_signal_dbm'] = snr_info.text
                if snr_info.tag == 'max_noise_dbm':
                    wn['max_noise_dbm'] = snr_info.text

        # Loop through GPS information
        if lev1.tag == 'gps-info':
            for gps_info in lev1:
                if gps_info.tag == 'peak-lat':
                    wn['peak_lat'] = gps_info.text
                if gps_info.tag == 'peak-lon':
                    wn['peak_lon'] = gps_info.text
    # select appropriate text for encryption field
    wn['encryption'] = populate_encryption(wn['placeholder_encryption'])

    print("Found infrastructure network with BSSID: %s - encryption: %s" \
             % (wn['bssid'], wn['encryption']))
    return wn

# Create an empty network dictionary with all needed keys
def make_net_dict():
    keys = ['wn_num',
            'first_seen',
            'last_seen',
            'ssid_type',
            'max_speed',
            'packets',
            'beaconrate',
            'wps',
            'wps_manuf',
            'dev_name',
            'model_name',
            'model_num',
            'placeholder_encryption',
            'encryption',
            'ssid_wpa_version',
            'cloaked',
            'essid',
            'bssid',
            'manuf',
            'channel',
            'maxseenrate',
            'max_signal_dbm',
            'max_noise_dbm',
            'clients',
            'peak_lat',
            'peak_lon']
    net_dict = {key: None for key in keys}
    return net_dict

# based in the entries in placeholder_encryption, return correct text
def populate_encryption(placeholder_list):
    encryption = 'UNKNOWN'

    if 'WEP' in placeholder_list and 'WPA' in placeholder_list:
        encryption = 'WEP + WPA'
    elif 'WEP' in placeholder_list:
        encryption = 'WEP'
    elif 'WPA+TKIP' in placeholder_list and 'WPA+PSK' in placeholder_list:
        encryption = 'WPA+TKIP/PSK'
    elif 'WPA+TKIP' in placeholder_list:
        encryption = 'WPA+TKIP'
    elif 'WPA+PSK' in placeholder_list:
        encryption = 'WPA+PSK'
    elif 'WPA+AES-CCM' in placeholder_list:
        encryption = 'WPA-MGT'
    else:
        encryption = 'OPEN'

    return encryption

## Own Code Contributions 

# Function for creating translation dictionary between essid and bssid
def e2bSSID_netlist(netList):
    e2bSSID = {}
    # 1st level is each wiresless network seen
    for wirelessNetwork in netList:
        # Each entry is a dict
        if debugBit != 0:
            print("Wireless Network:\t{0}\n\tType:\t{1}".format(wirelessNetwork,type(wirelessNetwork)))
            print("\t\tESSID:\t\t{0}\n\t\t\tType:\t{1}".format(wirelessNetwork['essid'],type(wirelessNetwork['essid'])))
            print("\t\tESSID:\t\t{0}\n\t\t\tType:\t{1}".format(wirelessNetwork['essid'],type(wirelessNetwork['essid'])))
        e2bSSID[wirelessNetwork['essid']] = wirelessNetwork['bssid']
    return e2bSSID

# Function to return bssid for a given essid
def retBSSID(essid, e2bList):
    # Check if the essid is in the e2b list
    if essid in e2bList:
        return e2bList[essid]
    else:
        print("[-] ISSUE: Passed essid ({0}) not found in e2b list".format(essid))
        return "DE:AD:BE:EF:CA:FE"

# Function for creating heatmaps for ONLY a SPECIFIC AP(s; repeaters)
def draw_data_by_essid(ts, bssid, signal, lat, lon, walked_lon, walked_lat, essid, e2bList):
    # Create a grid of 1000x1000 for plotting the data (NOTE: This will affect output based on size and scale o GPS data & range)
    grid_x, grid_y = numpy.mgrid[min(walked_lon):max(walked_lon):1000j, min(walked_lat):max(walked_lat):1000j]

    # Only deal with the unique APs
    bssids = list(set(bssid))

    # Determine the BSSID that matches the passed ESSID
    targetBSSID = retBSSID(essid, e2bList)

    # Loop through each BSSID
    for s in bssids:
        # Set up arrays for capture
        points_lon = []
        points_lat = []
        values = []
        h = []

        # Apply all points on an intermediate surface so that the code can distinct points where they are, without reception
        for i in range(0, len(bssid)):
            if debugBit != 0:
                print("[!] Test:\n\tbssid[i]\t=\t{0}\n\ts\t=\t{1}\n\ttargBSSID\t=\t{2}".format(bssid[i], s, targetBSSID))
            if bssid[i] == s and bssid[i] == targetBSSID:   # Check that we are looking at the correct bssid entry AND it matches the essid value
                hc = hash((lon[i], lat[i]))
                if hc not in h:
                    points_lon.append(lon[i])
                    points_lat.append(lat[i])
                    values.append(float(signal[i]))
                    h.append(hash((lon[i], lat[i])))

        # Optional: apply -100dBm where data was not gatherd
        for i in range(0, len(walked_lon)):
            hc = hash((walked_lon[i], walked_lat[i]))
            if hc not in h:
                points_lon.append(lon[i])
                points_lat.append(lat[i])
                values.append(float(-100))
                h.append(hash((walked_lon[i], walked_lat[i])))

        # Interpolate the data
        grid = griddata((points_lon, points_lat), numpy.array(values), (grid_x, grid_y), method='cubic')

        if s == targetBSSID:    # Check that we are only saving the image of the desired ESSID/BSSID
            # Store the bitmap in the current folder
            plt.show()
            #plt.imsave('%s.png' % (essid), grid.T)      # NOTE: Want to save as the essid
            plt.imsave('%s.png' % (essid), grid)    # Change graph to no longer use the transpose of the collected data (causing shift of data?)

# Function for creating heatmap over map image
def draw_data_path(lat, lon, filename):
##
# NOTE: Other python code for the overlay of a mapping ontop of a geolocaiton image:
#           1) Create Scatter plot that contains GPS locations as part of a large scatter plot (NOTE: THIS WILL CREATE WALKED PATH)
#           2) Import an image of the map to be used for the overlay
#               -> In example pulls in a specifically loaded image
#               - TODO: Have code grab a google earth (or other map) image based on GPS coordinates from collected data
#           3) Bound the box of the map used to match the GPS coordinate range of the information collected
# FROM https://www.sparkfun.com/tutorials/403   -   the gpsmap.py code
##
    # (Loop?) for creating a scatter plot made of all the collected GPS coordinates
    #plt.scatter(x=lon, y=lat, c='r')        # NOTE: Error may be coming from trying to set 'x' and 'y' by default
                                            #   -> Ex: plt.plot(x=lon, y=lat) RETURNS an empty array    | Try graphing two on same graph
    plt.plot(lat, lon, 'bo')

    # Set the labels for the scatter plot
    plt.ylabel('Longitude')
    plt.xlabel('Latitude')
    plt.title('Walked Path')

    # Import the map image to use for the overlay [ATTEMPT FIRST WITH HEATMAP IMAGE]
    #im = plt.imread('ATT2WU4IhY.png')   # Grabbed from already produced data
    im = plt.imread(filename)

    # Rotate the image file
    #im = np.rot90(im, k=-1) # NOTE: -1 was in the wrong direction
    im = np.rot90(im, k=1)      # NOTE: Seeing some weird TRANSPOSE rotation happeneing, but ONLY in ONE set of data.... WTF?

    # Get the max and min values from the GPS data
    min_lat=min(lat)
    max_lat=max(lat)
    min_lon=min(lon)
    max_lon=max(lon)
    print("Collected mins and maxes:\n\tMin Lat:\t\t{0}\n\tMax Lat:\t\t{1}\n\tMin Lon:\t\t{2}\n\tMax Lon:\t\t{3}".format(min_lat, max_lat, min_lon, max_lon))

    # Bound the box of the map used for the overlay
    #implot = plt.imshow(im, extent=[max_lat, min_lat, min_lon, max_lon])    # NOTE: Changed min_lat and max_lat to fix mismatch between
                                                                            #       heatmap orientation and graphed data points
    #implot = plt.imshow(im, aspect='auto', extent=[min_lat, max_lat, min_lon, max_lon])
    #implot = plt.imshow(im, aspect='auto', extent=[min_lat, max_lat, max_lon, min_lon])     # WORKS!!! Atleast the image and path match (ONLY FOR HOUSE DATA)

    implot = plt.imshow(im, aspect='auto', origin='lower', extent=[min_lat, max_lat, min_lon, max_lon]) # Seems to work, but not for everything

    # Show the produced graph
    plt.show()      # WORKS!!! For showing the path one walked over the heatmap that is being examined



# Function for printing out key to heatmap colors

# Function for combining all heatmap infomration(????)  NOTE: Not sure that this is really needed for what trying to do

# Function to graph both the heatmap information (2D) and the path data
#   -> NOTE: Needs to be able to hanlde multiple files of the same kind
def draw_path_data_by_essid(ts, bssid, signal, lat, lon, walked_lon, walked_lat, essid, e2bList):
    # Create a grid of 1000x1000 for plotting the data (NOTE: This will affect output based on size and scale o GPS data & range)
    grid_x, grid_y = numpy.mgrid[min(walked_lon):max(walked_lon):1000j, min(walked_lat):max(walked_lat):1000j]      # NOTE: Need to adjust so that mapping is NOT
                                                                                                                    #           forced into 1000x1000 grid

    # Only deal with the unique APs
    bssids = list(set(bssid))

    # Determine the BSSID that matches the passed ESSID
    targetBSSID = retBSSID(essid, e2bList)

    # Variables for storing path data
    x_lat = []
    y_lon = []

    # Loop through each BSSID
    for s in bssids:
        # Set up arrays for capture
        points_lon = []
        points_lat = []
        values = []
        h = []

        # Apply all points on an intermediate surface so that the code can distinct points where they are, without reception
        for i in range(0, len(bssid)):
            if debugBit != 0:
                print("[!] Test:\n\tbssid[i]\t=\t{0}\n\ts\t=\t{1}\n\ttargBSSID\t=\t{2}".format(bssid[i], s, targetBSSID))
            if bssid[i] == s and bssid[i] == targetBSSID:   # Check that we are looking at the correct bssid entry AND it matches the essid value
                hc = hash((lon[i], lat[i]))
                if hc not in h:
                    points_lon.append(lon[i])
                    points_lat.append(lat[i])
                    values.append(float(signal[i]))
                    h.append(hash((lon[i], lat[i])))
                    #plt.scatter(x=lat[i], y=lon[i], c='black')     # NOTE: Causes the scatter plot to become colored dots
                    x_lat.append(lat[i])
                    y_lon.append(lon[i])

        # Optional: apply -100dBm where data was not gatherd
        for i in range(0, len(walked_lon)):
            hc = hash((walked_lon[i], walked_lat[i]))
            if hc not in h:
                points_lon.append(lon[i])
                points_lat.append(lat[i])
                values.append(float(-100))
                h.append(hash((walked_lon[i], walked_lat[i])))
                #plt.scatter(x=lat[i], y=lon[i], c='black')        # NOTE: Still causes plot to become colored dots instead of a full image
                x_lat.append(lat[i])
                y_lon.append(lon[i])

        # Interpolate the data
        grid = griddata((points_lon, points_lat), numpy.array(values), (grid_x, grid_y), method='cubic')

        if s == targetBSSID:    # Check that we are only saving the image of the desired ESSID/BSSID
            # Attempt to rotate the grid data
            grid = np.rot90(grid, k=-1)
            # Store the bitmap in the current folder
            plt.show()
            #plt.imsave('%s.png' % (essid), grid.T)      # NOTE: Want to save as the essid
            plt.imsave('%s.png' % (essid), grid)    # Change graph to no longer use the transpose of the collected data (causing shift of data?)

            ##
            # Shifted in since it was grpahing for each BSSID seen
            ##
            # Attempt to graph all data
            plt.subplot(2, 1, 1)
            plt.plot(x_lat, y_lon, 'r.-')
            plt.xlabel('Latitude')
            plt.ylabel('Longitude')
            plt.title('Tale of Two Subplots')
    
            plt.subplot(2, 1, 2)
            plt.xlabel('Grid X Index')
            plt.ylabel('Grid Y Index')
            plt.title('Heatmap Grid')
            plt.imshow(grid)

            plt.show()
    print("[!] Test:\n\tbssid[i]\t=\t{0}\n\ts\t\t=\t{1}\n\ttargBSSID\t=\t{2}\n\n\tMin Lat:\t{3}\n\tMax Lat:\t{4}\n\tMin Lon:\t{5}\n\tMax Lon:\t{6}".format(bssid[i], s, targetBSSID, min(walked_lat), max(walked_lat), min(walked_lon), max(walked_lon)))

# Function for calling code sequence to produce heatmap + path [OLD WAY]
def path_map_combine(essidName, netFile, gpsFile):
    netList, clientList = load_nets_from_xml(netFile)
    e2bList = e2bSSID_netlist(netList)
    wifiName = essidName
    fileName = essidName + '.png'
    #draw_data_by_essid(*parse_xml(gpsFile), wifiName, e2bList)
    #ts, bssid, signal, lat, lon, walked_lon, walked_lat = parse_xml(gpsFile)
    #draw_data_path(lat, lon, fileName)
    draw_path_data_by_essid(*parse_xml(gpsFile), wifiName, e2bList)
    ts, bssid, signal, lat, lon, walked_lon, walked_lat = parse_xml(gpsFile)
    draw_data_path(lat, lon, fileName)

# From Github code:

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Usage %s << /path/to/Kismet.gpsxml >>" % (sys.argv[0]))
		sys.exit(-1)
	
	draw_data(*parse_xml(sys.argv[1]))

'''
 Example commands using code to produce wifi heatmap with overlapped GPS path graph

 1) netList, clientList = load_nets_from_xml('./Kismet-20181218-17-04-16-1.netxml')
    - Loads in the XML file to produce:
        a) netList
        b) clientList
 2) e2bList = e2bSSID_netlist(netList)
    - Creates a list that translates between ESSID and BSSID values
 3) wifiName = 'ATT2WU4IhY_vap'
    - Sets the name (ESSID) 
 4) draw_data_by_essid(*parse_xml('./Kismet-20181218-17-04-16-1.gpsxml'), wifiName, e2bList)
    - Creates an image (e.g. 'wifiName'.png) that is the heatmap for a specific ESSID (e.g. wifiName)
 5) ts, bssid, signal, lat, lon, walked_lon, walked_lat = parse_xml('./Kismet-20181218-17-04-16-1.gpsxml')
    - Creates arrays (ts, bssid, signal, lat, lon, walked_lon, walked_lat) for use with later code
 6) draw_data_path(lat, lon, filename)
    - Super imposes the GPS path data over the previously produced wifi heatmap.  NOTE: 'filename' should be the same as the file produced
        by 'draw_data_by_essid()'

 Fixed Run:
 plt.plot(lat, lon, 'bo')
 plt.xlabel('Longitude')
 plt.ylabel('Latitude')
 plt.title('Test Path')
 im = plt.imread(fileName)
 implot = plt.imshow(im, extent=[min_lat, max_lat, min_lon, max_lon])
 plt.show()

 House Test:
 path_map_combine('BEC450','./Kismet-20181228-15-08-56-1.netxml','Kismet-20181228-15-08-56-1.gpsxml')
'''
