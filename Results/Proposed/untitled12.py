fil=open("output.txt","r")
gil=open("intersect.txt","w")
try:
    while True:
        string=fil.readline()
        while (string[0]!="/"):
            gil.write(string)
            string=fil.readline()
        for i in range(5):
            point=open(f"iterate{i}.txt","w")
            fil.readline()
            string=fil.readline()
            while (string[0]!="/"):
                point.write(string)
                string=fil.readline()
            point.close()
            
except:
    fil.close()