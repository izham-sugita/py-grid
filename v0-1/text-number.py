


txt = "50"
a = "40"
a = a.zfill(3)
print(a)
x = txt.zfill(4)
print(x)
y = x.lstrip("0")
print(y)

filename = "naca0012-101-051-030-070-input.csv"
print(filename)
imax = filename[9:12]
imax = imax.lstrip("0")
imax = int(imax)
print( imax, type(imax))

jmax = filename[13:16]
jmax = jmax.lstrip("0")
jmax = int(jmax)
print( jmax, type(jmax))

tail1 = filename[17:20]
tail2 = filename[21:24]
print(tail1, tail2)

tail1 = tail1.lstrip("0")
tail2 = tail2.lstrip("0")

tail1 = int(tail1)
tail2 = int(tail2)
print(tail1, tail2, type(tail1), type(tail2) )
