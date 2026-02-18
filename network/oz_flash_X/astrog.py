#! /usr/bin/env python3


Initial = True



#with open('./in/rate-v3.z897c') as f:
#    for line in f:


with open('./in/rate-v3.z897c') as f:
    lines_in = [ line.rstrip("\n") for line in f ]




i_line = 0
reac_dat = []
header = []

line0  = ''
num_reac = 0
Initial = True

while( i_line < (len(lines_in) - 1) ):

    line = lines_in[i_line]

    dat = line.split()

    #print(i_line)

    if   dat[0]  == 'data':
        ### new reaction type
        if not Initial:
            header[-1][1] = num_reac

        header.append([dat[2], 0])

        Initial = False

        num_reac = 0
        i_line += 1
    elif dat[-1] == 'iaa':
        ### iaa (brussel data table)
        #print('iaa rate')



        table =      lines_in[i_line + 2] + ' ' + lines_in[i_line + 3] \
            + ' '  + lines_in[i_line + 4] + ' ' + lines_in[i_line + 5]

        reac_dat.append(['iaa', lines_in[i_line], lines_in[i_line + 1], \
                         qv, 0, table])

        num_reac += 1

        i_line += 6
    elif line == line0:
        ### resonance reaction fitting
        #print('resonance')


        a_fit =     lines_in[i_line + 2].split()[0] \
            + ' ' + lines_in[i_line + 2].split()[1] \
            + ' ' + lines_in[i_line + 2].split()[2] \
            + ' ' + lines_in[i_line + 2].split()[3] \
            + ' ' + lines_in[i_line + 3].split()[0] \
            + ' ' + lines_in[i_line + 3].split()[1] \
            + ' ' + lines_in[i_line + 3].split()[2]

        #print(reac_dat[-1])
        reac_dat[-1][4] += 1
        reac_dat[-1][5].append(a_fit)

        #print(reac_dat[-1])


        i_line += 4
    else:
        line0 = line

        a_fit =     lines_in[i_line + 2].split()[0] \
            + ' ' + lines_in[i_line + 2].split()[1] \
            + ' ' + lines_in[i_line + 2].split()[2] \
            + ' ' + lines_in[i_line + 2].split()[3] \
            + ' ' + lines_in[i_line + 3].split()[0] \
            + ' ' + lines_in[i_line + 3].split()[1] \
            + ' ' + lines_in[i_line + 3].split()[2]

        qv      = lines_in[i_line + 3].split()[3]

        reac_dat.append(['reaclib', lines_in[i_line], lines_in[i_line + 1], \
                         qv, 0, [a_fit]])

        num_reac += 1
        i_line += 4


header[-1][1] = num_reac


###check

if False:

    for j in range(len(reac_dat)):
        if reac_dat[j][4] > 0 and False:
            print(reac_dat[j])

        if reac_dat[j][0] == 'iaa' and True:
            print(reac_dat[j])




### output

Out = open('./out/ntwk.dat', 'w')

ireac = - 1

for itype in range(len(header)):
    Out.write(' data of  ' + header[itype][0] + ' reactions\n')


    for ireac_in in range(header[itype][1]):
        #print(reac_dat[ireac][0])

        ireac += 1

        if   reac_dat[ireac][0] == 'iaa':
            Out.write(reac_dat[ireac][1] + '\n')
            Out.write(reac_dat[ireac][2] + '\n')

            table = reac_dat[ireac][5]
            table = table.split()
            for j in range(len(table)):
                Out.write('{:>11}'.format(table[j]))
                if j % 8 == 7:
                    Out.write('\n')
            Out.write('\n')

        elif reac_dat[ireac][0] == 'reaclib':
            #print(reac_dat[ireac])
            Out.write(reac_dat[ireac][1] + '\n')
            Out.write(reac_dat[ireac][2] + '\n')
            qv = float(reac_dat[ireac][3])
            for j in range(1 + reac_dat[ireac][4]):
                a_fit = reac_dat[ireac][5][j]
                a_fit = a_fit.split(' ')
                Out.write('{0:13.5e}'.format(float(a_fit[0])))
                Out.write('{0:13.5e}'.format(float(a_fit[1])))
                Out.write('{0:13.5e}'.format(float(a_fit[2])))
                Out.write('{0:13.5e}\n'.format(float(a_fit[3])))
                Out.write('{0:13.5e}'.format(float(a_fit[4])))
                Out.write('{0:13.5e}'.format(float(a_fit[5])))
                Out.write('{0:13.5e}'.format(float(a_fit[6])))
                Out.write('{0:13.5e}\n'.format(qv))

        else:
            print('error #1: please check')
            exit()




Out.write('  end   of data set                                \n\n')









exit()
