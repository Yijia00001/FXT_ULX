#!/bin/bash
## run specextract to generate spectrum 

cat obsid_list.txt | while read line;
do
    if [ `grep -c "$line" finish.txt` -ne '0' ]; then
        continue
    fi
    
    cd data/$line
    pwd
    inst_list=("a" "b")
    for inst in ${inst_list[@]};
    do
        cd ${inst}
        evt_file=`ls ../fxt/products/*_${inst}_*`
        # img_file=$inst.img

        source_num=`ls x*.reg | wc -l`
        source_index=1
        while [ $source_index -le ${source_num} ]
        do
            mkdir x${source_index}
            cd x${source_index}

            evt_file=`ls ../../fxt/products/*_${inst}_*`
            echo "

read events ${evt_file}
.
yes
filter region ../x${source_index}.reg
extract spectrum
save spectrum x${source_index}.pi clobber=yes
clear region
filter region ../bkg.reg
extract spectrum
save spectrum bkg.pi clobber=yes
exit
no" > xselect_para.xco

            xselect @xselect_para.xco
            fxtarfgen specfile=x${source_index}.pi expfile=../expo_${inst}-without_vign.img outfile=x${source_index}.arf clobber=yes
            fxtrmfgen specfile=x${source_index}.pi outfile=x${source_index}.rmf clobber=yes
            # grppha infile=x${source_index}.pi outfile=x${source_index}_grp.pi chatter=0 \
            # comm=" group min 15 & chkey RESPFILE x${source_index}.rmf & chkey ANCRFILE x${source_index}.arf \
            # & chkey BACKFILE bkg.pi & exit" clobber=yes
            #break 1
            cd ..
            let source_index++
        done
        #break 1
        cd ..
    done
    cd ../..
    #break 1
done
