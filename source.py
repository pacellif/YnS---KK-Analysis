import paramiko
import getpass
import filereader as fr
import ROOT

# SSH connection parameters
hostname = 'lxplus8.cern.ch'
port = 22  # default SSH port
username = 'fpacelli'
password = getpass.getpass("Password: ")  # or use key-based authentication

# Path to the file on the remote server
remote_file_path = '/afs/cern.ch/user/f/fpacelli/private/ListY2SPhiRun2.txt'

# Create an SSH client
ssh = paramiko.SSHClient()

# Automatically add the server's host key (this is insecure; see the Paramiko documentation for how to do it securely)
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

# Connect to the server
ssh.connect(hostname, port, username, password)

# Open the file using SFTP (SSH File Transfer Protocol)
with ssh.open_sftp() as sftp:
    # Open the remote file
    remote_file = sftp.file(remote_file_path, 'r')

    # Read the contents of the file
    file_contents = remote_file.readlines()

    # Do something with the file contents (e.g., print or process)
    
for i in range(len(file_contents)):			
	file_contents[i] = file_contents[i].replace("\n", "")

sample = file_contents[:4]    

dataY2SKK = ROOT.RDataFrame("rootuple/CandidateTree", sample)
if not dataY2SKK:
	print("Error....")
	exit()
dataY2S = ROOT.RDataFrame("rootuple/UpsTree", sample)
print("fin qui tuto bene?")

# Close the SSH connection
#ssh.close()





















