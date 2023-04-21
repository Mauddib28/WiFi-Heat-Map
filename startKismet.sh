#!/bin/sh

###
# The purpose of this script is to begin running kismet in a custom directory 
###

## Variables
debugBit=1

## Begin script
echo -e "[*] Preparing to kismet"

# Ask the user for directory name to save the data under
read -p "What directory would you like to store the data in?" userDir
if [[ debugBit -eq 1 ]]; then 
	echo -e "\tYou want the directory: $userDir"
fi

# Check to see if the directory does not exist
if [ ! -d $userDir ]; then		# Does the directory not exist?
	if [[ debugBit -eq 1 ]]; then
		echo -e "[-] User supplied directory does not exist.... Creating it"
	mkdir $userDir
fi

# Move into the user chosen directory
cd $userDir

# Begin running kismet within a tmux terminal
#tmux new-session 'kismet' \;

# Begin the kismet_server as a background processs
#	NOTE: Need to be sure that running as root!!
sudo kismet_server --daemonize		# Seems that data collection starts with this

# Will eventually need a way to kill the background process....
#	-> Hopefully better than 'kill -9 kismet_server'
