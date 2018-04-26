
#               (25/04/2018)
#            @Joanna Suliga Jan 2018
#                    after
#            @Joy Bhattacharjee Jun 2017
#

#-------------------------------------------------------
# - - - MODULES AND WORKING DIRECTORY - - - - - - - - -
#-------------------------------------------------------

#importing all the modules
from ftplib import FTP
import os,os.path


directory = os.path.dirname(os.path.realpath(__file__))
dir_input_maps = directory + "\\probaV_download"

# ftp server default address
os.chdir(dir_input_maps)
ftp = FTP('ftp.vito-eodata.be')

#-------------------------------------------------------
# - - - USER INPUT - - - - - - - - - - - - - - - - - - -
#-------------------------------------------------------
# provide username and password
username=raw_input("Provide your username: ")
pswrd=raw_input("Provide your password: ")
print 'Logging in.'
ftp.login(username,pswrd)

# provide ftp directory
vito_download = raw_input('Provide order name given in email from VITO (ex.M0169702) : ')

# taking the list of the directory    
print 'Changing to ' + vito_download 
ftp.cwd(vito_download )
ftp.retrlines('LIST')
print 'Downloading files'

# get filenames within the folder
filenames = ftp.nlst() # get filenames within the directory

#Creating folder list
folder_list=[]
for i in filenames:
    all_list=ftp.nlst(i)
    folder_list.append(all_list)

# saving the files in local directory
for name in folder_list:
    for i in name:
        local_filename = os.path.join(dir_input_maps,i[34:])
        file = open(local_filename, 'wb')
        ftp.retrbinary('RETR '+ i, file.write)
        file.close()
ftp.quit() 

print "Done!"
