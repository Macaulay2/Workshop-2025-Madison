#!/bin/bash
set -e  # Stop on error

echo "========== STARTING SETUP =========="

# Making sure Ubuntu packages are up to date. 
echo "Updating package lists..."
sudo apt -y update

# Installing dependencies necessary to install Macaulay2
echo "Installing dependencies..."
sudo apt install -y software-properties-common curl

# Setting default options fro the postfix config pop-up so we don't need manual input. 
echo "postfix postfix/main_mailer_type select No configuration" | sudo debconf-set-selections
echo "postfix postfix/mailname string localhost" | sudo debconf-set-selections

# Add the repository PPA for Macaulay2
sudo add-apt-repository -y ppa:macaulay2/macaulay2
# Updates the repository to make sure that it's up to date
sudo apt -y update
# Installing macaulay2, making sure that any frontend interaction is suppressed with the previously specificed preset configs. 
sudo DEBIAN_FRONTEND=noninteractive apt install -y macaulay2 

# Adding aliases so that it's easier to call m2 via the terminal, whether bash or zshrc 
echo "alias m2='M2'" >> ~/.bashrc
echo "alias m2='M2'" >> ~/.zshrc