#! /usr/bin/env python3


Initial = True


header = []
line0  = ''

#with open('./in/rate-v3.z897c') as f:
#    for line in f:


with open('./in/rate-v3.z897c') as f:
    lines_in = f.readlines()


i_line = 0
reac_dat = []

while( i_line < (len(lines_in) - 1) ):

    line = lines_in[i_line]

    dat = line.split()

    print(i_line)

    if   dat[0]  == 'data':
        ### new reaction type
        header.append(dat[1])
        i_line += 1
    elif dat[-1] == 'iaa':
        ### iaa (brussel data table)
        print('iaa rate')

        i_line += 6
    elif line == line0:
        ### resonance reaction fitting
        print('resonance')


        a_fit =     lines_in[i_line + 2].split()[0] + ' ' + lines_in[i_line + 2].split()[1] \
            + ' ' + lines_in[i_line + 2].split()[2] + ' ' + lines_in[i_line + 2].split()[3] \
            + ' ' + lines_in[i_line + 3].split()[0] + ' ' + lines_in[i_line + 3].split()[1] \
            + ' ' + lines_in[i_line + 3].split()[2]

        print(reac_dat[-1])

        reac_dat[-1][4] += 1
        reac_dat[-1][5].append(a_fit)

        print(reac_dat[-1])

        exit()


        i_line += 4
    else:
        line0 = line

        a_fit =     lines_in[i_line + 2].split()[0] + ' ' + lines_in[i_line + 2].split()[1] \
            + ' ' + lines_in[i_line + 2].split()[2] + ' ' + lines_in[i_line + 2].split()[3] \
            + ' ' + lines_in[i_line + 3].split()[0] + ' ' + lines_in[i_line + 3].split()[1] \
            + ' ' + lines_in[i_line + 3].split()[2]

        qv      = lines_in[i_line + 3].split()[3]

        reac_dat.append(['reaclib', lines_in[i_line], lines_in[i_line + 1], qv, 0, [a_fit]])
        i_line += 4

print(reac_dat[1010])

exit()
