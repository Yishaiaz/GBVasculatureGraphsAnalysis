import os

try:
    with open('job_conf_file.txt', 'r') as f:
        print("FOUND FILE!, the content is:")
        print(f.readlines())
except FileNotFoundError as e:
    print("could not find file, but here is the cwd listdir")
    print(f"CWD:{os.getcwd()}")
    print(os.listdir(os.getcwd()))
