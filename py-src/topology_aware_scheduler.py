import struct

# Sample list of integers
int_list = [42, 13, -7, 0, 66545535]

f = open('integers.bin', 'wb')

for ii in int_list:

    f.write((ii).to_bytes(8, byteorder='little', signed=True))

f.close()




f = open('integers.bin', 'rb')



int_list = []
while True: 
    chunk = f.read(8)
    if not chunk: break 

    int_list.append(struct.unpack('q', chunk)[0])


f.close()

print(int_list)

