import os 

path = os.getcwd() 
print(path)
new_path = os.path.join(path,'pictures') + os.sep
print(new_path) 