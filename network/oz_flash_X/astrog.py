#! /usr/bin/env python3



Ra26     = True
Ra26_mod = True

Ra26_all = True
#Ra26_all = False #only rath & ths8

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

Out.close()

print('')


#####  Rauscher 2026 ###################################3
if Ra26:

    print(' ----- read Rauscher 2026 data -----')

    ### p-induced
    ra26_pg = []; ra26_pn = []; ra26_pa = []
    with open('../../lib_data/ra26/fits_p.dat') as f:
        i_cnt = 0
        for line in f:
            line = line.rstrip('\n')
            i_cnt += 1
            if i_cnt % 3 == 1:
                dat = line.split()
                nuc1 = dat[0]
                nuc2 = dat[3]
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
                    ra26_pg.append([nuc1, nuc2, line])
                elif rtype == 'pn':
                    ra26_pn.append([nuc1, nuc2, line])
                elif rtype == 'pa':
                    ra26_pa.append([nuc1, nuc2, line])
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
            line = line.rstrip('\n')
            i_cnt += 1
            if i_cnt % 3 == 1:
                dat = line.split()
                nuc1 = dat[0]
                nuc2 = dat[3]
                #print(dat[1],dat[2])
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
                    ra26_ng.append([nuc1, nuc2, line])
                elif rtype == 'np':
                    ra26_np.append([nuc1, nuc2, line])
                elif rtype == 'na':
                    ra26_na.append([nuc1, nuc2, line])
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
            line = line.rstrip('\n')
            if i_cnt % 3 == 1:
                dat = line.split()
                nuc1 = dat[0]
                nuc2 = dat[3]
                #print(dat[1],dat[2])
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
                    ra26_ag.append([nuc1, nuc2, line])
                elif rtype == 'an':
                    ra26_an.append([nuc1, nuc2, line])
                elif rtype == 'ap':
                    ra26_ap.append([nuc1, nuc2, line])
            else:
                if   rtype == 'ag':
                    ra26_ag[-1].append(line)
                elif rtype == 'an':
                    ra26_an[-1].append(line)
                elif rtype == 'ap':
                    ra26_ap[-1].append(line)

    print(' ---------- finished ----------' + '\n')



####################################################

print(' ----- read rate-v3.z897c data -----')

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

print(' ---------- finished ----------' + '\n')




if Ra26 and Ra26_mod:
    ### find reaction rates included in Ra26

    ireac = -1  ### counter for reac_dat[]

    for itype in range(len(header)):

        #Reverse = False
        Mod = True
        if   header[itype][0] == 'pg':
            rt_for = ra26_pg.copy()
            rt_rev = []
            i_nuc  = 1
        if   header[itype][0] == 'gp':
            rt_for = []
            rt_rev = ra26_pg.copy()
            i_nuc  = 0
        elif header[itype][0] == 'pn':
            rt_for = ra26_pn.copy()
            rt_rev = []
            i_nuc  = 1
        elif header[itype][0] == 'np':
            rt_for = ra26_np.copy()
            rt_rev = []
            i_nuc  = 1
        elif header[itype][0] == 'pa':
            rt_for = ra26_pa.copy()
            rt_rev = []
            i_nuc  = 1
        elif header[itype][0] == 'ap':
            rt_for = ra26_ap.copy()
            rt_rev = []
            i_nuc  = 1
        elif header[itype][0] == 'ng':
            rt_for = ra26_ng.copy()
            rt_rev = []
            i_nuc  = 1
        elif header[itype][0] == 'gn':
            rt_for = []
            rt_rev = ra26_ng.copy()
            i_nuc  = 0
        elif header[itype][0] == 'ag':
            rt_for = ra26_ag.copy()
            rt_rev = []
            i_nuc  = 1
        elif header[itype][0] == 'ga':
            rt_for = []
            rt_rev = ra26_ag.copy()
            i_nuc  = 0
        elif header[itype][0] == 'na':
            rt_for = ra26_na.copy()
            rt_rev = []
            i_nuc  = 1
        elif header[itype][0] == 'an':
            rt_for = ra26_an.copy()
            rt_rev = []
            i_nuc  = 1
        else:
            ireac += header[itype][1]
            Mod = False
            print(' - no modify', header[itype][0])


        if Mod:

            print('     - modifying: ' + header[itype][0])

            for ireac_in in range(header[itype][1]):

                ireac += 1

                if   reac_dat[ireac][0] == 'iaa':
                    #print('iaa')
                    continue
                elif reac_dat[ireac][0] == 'reaclib':

                    if   i_nuc == 0:
                        nuc0  = reac_dat[ireac][1][0:5].strip()
                    elif i_nuc == 1:
                        nuc0  = reac_dat[ireac][1][5:10].strip()
                    else:
                        exit('error: i_nuc')

                    nuc0  = nuc0.lower()
                    label = reac_dat[ireac][1][41:45].strip()

                    Reverse = False
                    if reac_dat[ireac][1][39:40] == 'v':
                        Reverse = True

                    Mod_rate = False
                    if   Ra26_all:
                        Mod_rate = True
                    elif label == 'rath' or label == 'ths8':
                        Mod_rate = True

                    if Mod_rate:
                        #print(label)
                        for j in range(len(rt_for)):
                            if nuc0 == rt_for[j][0]:
                                reac_dat[ireac][1]\
                                    = reac_dat[ireac][1][0:41] + 'ra26'

                                if reac_dat[ireac][4] != 0:
                                    reac_dat[ireac][4] = 0
                                    reac_dat[ireac][-1] = 'tmp'

                                if Reverse:
                                    if len(rt_rev) > 1:
                                        a_tmp = rt_rev[j][3].lstrip()
                                    else:
                                        print('####### problem', reac_dat[ireac])
                                else:
                                    a_tmp = rt_for[j][2].lstrip()

                                a_tmp = ' '.join(a_tmp.split())

                                reac_dat[ireac][-1] = [a_tmp]

                else:
                    print('error #2')
                    exit()



#print(reac_dat[ireac])


#exit('in prep.')


###### output ##################################################

print(' ----- out new rate data -----')

Out = open('./out/rate.dat', 'w')

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
            qv = float(reac_dat[ireac][3])
            for j in range(1 + reac_dat[ireac][4]):
                a_fit = reac_dat[ireac][5][j]
                a_fit = a_fit.split(' ')

                #print(a_fit)
                Out.write(reac_dat[ireac][1] + '\n')
                Out.write(reac_dat[ireac][2] + '\n')

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




Out.write('  end   of data set\n')


print(' ---------- finished ----------' + '\n')






exit()
