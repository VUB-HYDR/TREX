#           ProbaV - downloading
#                 (30/05/2018)
#-------------------------------------------------------
# - - - MODULES AND WORKING DIRECTORIES - - - - - - - - -
#-------------------------------------------------------

from ftplib import FTP
import os,os.path

directory = os.path.dirname(os.path.realpath(__file__))
dir_input_maps = directory + "\\probaV_download"
os.chdir(dir_input_maps)
ftp = FTP('ftp.vito-eodata.be')

#-------------------------------------------------------
# - - - USER INPUT - - - - - - - - - - - - - - - - - - -
#-------------------------------------------------------
# For more information regarding downloading images please explore ReadMe file.
# Script will ask to provide username and password used at www.vito-eodata.be
username=raw_input("Provide your username: ")
pswrd=raw_input("Provide your password: ")
print 'Logging in.'
ftp.login(username,pswrd)

# and provide ftp directory.
vito_download = raw_input('Provide order name given in the email from VITO (mailing@vito.be) for example M0169702) : ')

#-------------------------------------------------------
# - - - DOWNLOADING - - - - - - - - - - - - - - - - - - -
#-------------------------------------------------------
ftp.cwd(vito_download )
ftp.retrlines('LIST')
print 'Downloading files'

filenames = ftp.nlst()
folder_list=[]
for i in filenames:
    all_list=ftp.nlst(i)
    folder_list.append(all_list)

for name in folder_list:
    for i in name:
        local_filename = os.path.join(dir_input_maps,i[34:])
        file = open(local_filename, 'wb')
        ftp.retrbinary('RETR '+ i, file.write)
        file.close()
ftp.quit() 

print "Downloading complete."
