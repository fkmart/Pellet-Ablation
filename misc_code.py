#This file is a congregation of various pieces of code I regularly

#Code set a filepath, makes a list of files whose title starts with 'title' and deletes them
path = os.path.join(os.path.expanduser('~'), 'Pictures')
filelist = [ f for f in os.listdir(path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(path, f))