
import os, sys

apui = []
syui = []
noui = []

root = sys.argv[1]
for dirName, subdirList, fileList in os.walk(root):
    if True in [dirName.find(s)>=0 for s in ("__","build",".git",".vscode") ]:
        continue
    for f in fileList:
        if f.endswith(".py"):
            p = os.path.join(dirName, f)
            s = open( p , "r" ).read()
            if "optparse" in s or "add_option" in s:
                apui.append( p )
            elif "argv" in s:
                syui.append( p )
            else:
                noui.append( p )

print("Neither:")
for f in noui:
    print(f)


print("argv:")
for f in syui:
    print(f)    

print("optparse:")
for f in apui:
    print(f)
