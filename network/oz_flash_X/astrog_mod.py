#! /usr/bin/env python3

import nuc_data as nd


Rate_mod = True # switch ON/OFF rate modification;
                # Original rates are outputted if Rate_mod = False.


#### Which reaction rates are updated by Ra26 (Rauscher 2026)?
Ra26     = True # must be true

Ra26_all = False # all the rates in Ra26 (not only theory rates)

n_selec = 4 # additional slection rules:
            #   0: all
            #   1: SnSbTe  test
            #   2: only pg (& gp)
            #   3: only alpha-induced reactions
            #   4: only n-induced reactions


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




if Ra26 and Rate_mod:
    ### find reaction rates included in Ra26

    # Each Ra26 entry has the structure
    #
    #     [target_nucleus, product_nucleus, forward_fit, reverse_fit]
    #
    # For particle-exchange reactions, the desired rate may be stored either
    # directly or as the reverse rate of a reaction in another input file.
    #
    # Each tuple below is:
    #
    #     (rate_list, nucleus_index, fit_index, match_type)
    #
    # A direct match compares nuc0 with entry[0] and uses entry[2].
    # A reverse match compares nuc0 with entry[1] and uses entry[3].

    reaction_search = {
        # Radiative capture and photodisintegration
        'pg': [
            (ra26_pg, 0, 2, 'direct'),
        ],
        'gp': [
            (ra26_pg, 1, 3, 'reverse'),
        ],
        'ng': [
            (ra26_ng, 0, 2, 'direct'),
        ],
        'gn': [
            (ra26_ng, 1, 3, 'reverse'),
        ],
        'ag': [
            (ra26_ag, 0, 2, 'direct'),
        ],
        'ga': [
            (ra26_ag, 1, 3, 'reverse'),
        ],

        # Particle-exchange reactions.
        # Search the same direction first, then the reverse direction
        # stored in the corresponding reaction list from another file.
        'pn': [
            (ra26_pn, 0, 2, 'direct'),
            (ra26_np, 1, 3, 'reverse'),
        ],
        'np': [
            (ra26_np, 0, 2, 'direct'),
            (ra26_pn, 1, 3, 'reverse'),
        ],
        'pa': [
            (ra26_pa, 0, 2, 'direct'),
            (ra26_ap, 1, 3, 'reverse'),
        ],
        'ap': [
            (ra26_ap, 0, 2, 'direct'),
            (ra26_pa, 1, 3, 'reverse'),
        ],
        'na': [
            (ra26_na, 0, 2, 'direct'),
            (ra26_an, 1, 3, 'reverse'),
        ],
        'an': [
            (ra26_an, 0, 2, 'direct'),
            (ra26_na, 1, 3, 'reverse'),
        ],
    }

    # Position of the target nucleus in the first line of each rate entry.
    target_column = {
        'pg': 1,
        'gp': 0,
        'ng': 1,
        'gn': 0,
        'ag': 1,
        'ga': 0,
        'pn': 1,
        'np': 1,
        'pa': 1,
        'ap': 1,
        'na': 1,
        'an': 1,
    }

    ireac = -1  ### counter for reac_dat[]

    for itype in range(len(header)):

        reaction_type = header[itype][0]

        if reaction_type not in reaction_search:
            ireac += header[itype][1]
            print(' - no modify', reaction_type)
            continue

        print('     - modifying: ' + reaction_type)

        i_nuc = target_column[reaction_type]
        search_candidates = reaction_search[reaction_type]

        count_theory = 0
        count_direct = 0
        count_reverse = 0
        count_not_found = 0
        count_missing_fit = 0

        for ireac_in in range(header[itype][1]):

            ireac += 1

            if reac_dat[ireac][0] == 'iaa':
                continue

            if reac_dat[ireac][0] != 'reaclib':
                print('error #2')
                exit()

            if i_nuc == 0:
                nuc0 = reac_dat[ireac][1][0:5].strip()
            elif i_nuc == 1:
                nuc0 = reac_dat[ireac][1][5:10].strip()
            else:
                exit('error: i_nuc')

            nuc0 = nuc0.lower()
            label = reac_dat[ireac][1][41:45].strip()

            Mod_rate = False

            if Ra26_all:
                Mod_rate = True

            elif label in ('rath', 'ths8', 'thra'):

                if n_selec == 0:
                    Mod_rate = True

                elif n_selec == 1:
                    nuc_tmp = nd.iso_name(nuc0)
                    if nuc_tmp[1] < 47 and nuc_tmp[2] < 53:
                        Mod_rate = True

                elif n_selec == 2:
                    Mod_rate = reaction_type in ('pg', 'gp')

                elif n_selec == 3:
                    Mod_rate = reaction_type in (
                        'ap', 'pa', 'ag', 'ga'
                    )

                elif n_selec == 4:
                    Mod_rate = reaction_type in (
                        'ng', 'gn', 'np', 'pn'
                    )

                else:
                    print('under construction')
                    exit()

            if not Mod_rate:
                continue

            count_theory += 1

            found = False
            a_tmp = None
            match_type = None

            for rate_list, nuc_index, fit_index, search_type \
                    in search_candidates:

                for rate_entry in rate_list:

                    if nuc0 != rate_entry[nuc_index]:
                        continue

                    if len(rate_entry) <= fit_index:
                        count_missing_fit += 1
                        print(
                            '####### missing fit:',
                            reaction_type,
                            nuc0,
                            search_type
                        )
                        continue

                    a_tmp = rate_entry[fit_index].lstrip()
                    match_type = search_type
                    found = True
                    break

                if found:
                    break

            if not found:
                count_not_found += 1
                print(
                    '####### Ra26 rate not found:',
                    reaction_type,
                    nuc0,
                    label
                )
                continue

            if match_type == 'direct':
                count_direct += 1
            elif match_type == 'reverse':
                count_reverse += 1
            else:
                exit('error: unknown Ra26 match type')

            reac_dat[ireac][1] = \
                reac_dat[ireac][1][0:41] + 'ra26'

            # Remove additional resonance fits and retain only the Ra26 fit.
            reac_dat[ireac][4] = 0

            a_tmp = ' '.join(a_tmp.split())
            reac_dat[ireac][-1] = [a_tmp]

        print(
            '       theory={0}, direct={1}, reverse={2}, '
            'not found={3}, missing fit={4}'.format(
                count_theory,
                count_direct,
                count_reverse,
                count_not_found,
                count_missing_fit
            )
        )



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
