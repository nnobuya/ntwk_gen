#! /usr/bin/env python3



Ra26 = True

### part data

Out = open('./out/net_nz.dat', 'w')

nuc_ntwk = []; a_ntwk = []; n_ntwk = []; z_ntwk = []
with open('./in/part-v3.z897') as f:
    for line in f:
        dat = line.split()
        if len(dat) == 6:
            nuc_ntwk.append(dat[0])
            a_ntwk.append(int(dat[1]))
            n_ntwk.append(int(dat[2]))
            z_ntwk.append(int(dat[3]))

for j in range(len(nuc_ntwk)):
    Out.write('{0:5}{1:5}{2:5}\n'.format(z_ntwk[j], n_ntwk[j], a_ntwk[j]))



#####  Rauscher 2026 ###################################3
if Ra26:

    ### p-induced
    ra26_pg = []; ra26_pn = []; ra26_pa = []
    with open('../../lib_data/ra26/fits_p.dat') as f:
        i_cnt = 0
        for line in f:
            i_cnt += 1
            if i_cnt % 3 == 1:
                dat = line.split()
                nuc_tmp = dat[0]
                if   dat[2] == 'gam':
                    rtype = 'pg'
                elif dat[2] == 'n':
                    rtype = 'pn'
                elif dat[2] == 'he4':
                    rtype = 'pa'
                else:
                    print(dat)
                    exit('error: check fits_p')
            elif i_cnt % 3 == 2:
                if   rtype == 'pg':
                    ra26_pg.append([nuc_tmp, line])
                elif rtype == 'pn':
                    ra26_pn.append([nuc_tmp, line])
                elif rtype == 'pa':
                    ra26_pa.append([nuc_tmp, line])
            else:
                if   rtype == 'pg':
                    ra26_pg[-1].append(line)
                elif rtype == 'pn':
                    ra26_pn[-1].append(line)
                elif rtype == 'pa':
                    ra26_pa[-1].append(line)


    ### n-induced
    ra26_ng = []; ra26_np = []; ra26_na = []
    with open('../../lib_data/ra26/fits_n.dat') as f:
        i_cnt = 0
        for line in f:
            i_cnt += 1
            if i_cnt % 3 == 1:
                dat = line.split()
                nuc_tmp = dat[0]
                print(dat[1],dat[2])
                if   dat[2] == 'gam':
                    rtype = 'ng'
                elif dat[2] == 'p':
                    rtype = 'np'
                elif dat[2] == 'he4':
                    rtype = 'na'
                else:
                    print(dat)
                    exit('error: check fits_n')
            elif i_cnt % 3 == 2:
                if   rtype == 'ng':
                    ra26_ng.append([nuc_tmp, line])
                elif rtype == 'np':
                    ra26_np.append([nuc_tmp, line])
                elif rtype == 'na':
                    ra26_na.append([nuc_tmp, line])
            else:
                if   rtype == 'ng':
                    ra26_ng[-1].append(line)
                elif rtype == 'np':
                    ra26_np[-1].append(line)
                elif rtype == 'na':
                    ra26_na[-1].append(line)


    ### a-induced
    ra26_ag = []; ra26_an = []; ra26_ap = []
    with open('../../lib_data/ra26/fits_a.dat') as f:
        i_cnt = 0
        for line in f:
            i_cnt += 1
            if i_cnt % 3 == 1:
                dat = line.split()
                nuc_tmp = dat[0]
                print(dat[1],dat[2])
                if   dat[2] == 'gam':
                    rtype = 'ag'
                elif dat[2] == 'n':
                    rtype = 'an'
                elif dat[2] == 'p':
                    rtype = 'ap'
                else:
                    print(dat)
                    exit('error: check fits_a')
            elif i_cnt % 3 == 2:
                if   rtype == 'ag':
                    ra26_ag.append([nuc_tmp, line])
                elif rtype == 'an':
                    ra26_an.append([nuc_tmp, line])
                elif rtype == 'ap':
                    ra26_ap.append([nuc_tmp, line])
            else:
                if   rtype == 'ag':
                    ra26_ag[-1].append(line)
                elif rtype == 'an':
                    ra26_an[-1].append(line)
                elif rtype == 'ap':
                    ra26_ap[-1].append(line)





exit()






####

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
