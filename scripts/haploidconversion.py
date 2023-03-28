import sys
chromosome = int(sys.argv[1])
hap_to_dip = {'0':'0/0', '1':'1/1'}

with open(f'data/hartfield/atsweeps/analysis/byanez/VCF/Vcf_Chr{chromosome}_Neutral.vcf',"r") as file, open("/scratch/byanez/Vcf_Chr5_Neutral_Diploid.vcf","w") as s_file:
        for line in file:
                if line[0] == '#':
                        s_file.write(line)
                else:
                        rowstr = line.split('\t')[9:]
                        for i in range(len(rowstr)):
                                rowstr[i] = hap_to_dip[rowstr[i][0]] + rowstr[i][1:]
                                print(rowstr)
                        row = '\t'.join(str(line).split('\t')[:9]) + '\t' + '\t'.join(rowstr)
                        s_file.write(row)

