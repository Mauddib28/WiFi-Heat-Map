#!/bin/sh

###
# The purpose of this script is to check if the GPS module has a fix
#	NOTE:	This script assumes that the GPS module communicates through /dev/ttyAMA0
###

## Variables
debugBit=1
gpsDEVfile="/dev/ttyAMA0"

## Begin script
echo "[*] Checking for a GPS lock"

# Check if a GPS lock has been seen
#cat $gpsDEVfile | grep GPGGA | cut -d',' -f 7		# WORKS, but will be continuous... need something that is time restricted

# While loop to search for GPGGA entry from GPS module
echo "[*] Beginning loop to check for GPS fix"
while true; do
	if [ "$debugBit" -eq "1" ]; then
		echo "\t[*] About to read from $gpsDEVfile"
		echo "[*] Test reading from GPS file..."
	fi
	checkGPS=$(cat /dev/ttyAMA0 | grep -m 1 GPGGA)		# Not -c1, this counts the number of occurances....., use '-m 1'
	if [ "$debugBit" -eq "1" ]; then
		echo "\t[?] checkGPS value: $checkGPS"
	fi
	if [ "$(echo $checkGPS | cut -d',' -f7)" -ne "0" ]; then	# Check if the fix quality is non-zero
		echo "[+] GPS has a lock....."
		break
	fi
done

# Return information letting the user know that a GPS lock has been seen
echo "[*] GPS lock has been found"
exit 0
