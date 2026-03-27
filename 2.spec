#!/bin/bash
## run specextract to generate spectrum 



cat obsid_list.txt | while read line;
do
    if [ `grep -c "$line" finish.txt` -ne '0' ]; then
        continue
    fi
    cd data/$line

    mkdir bkg_spec
    cd bkg_spec
    evt_file=`ls ../fxt/products/*fxt_a*.fits`
    echo "

read events ${evt_file}
.
yes
filter region ../bkg.reg
extract spectrum
save spectrum bkg_a.pi clobber=yes
clear region
exit
no" > xselect_para_bkg.xco
    xselect @xselect_para_bkg.xco
    fxtarfgen specfile=bkg_a.pi expfile=../a/expo_a-without_vign.img outfile=bkg_a.arf clobber=yes
    fxtrmfgen specfile=bkg_a.pi outfile=bkg_a.rmf clobber=yes
    cd ..

    inst_list=("a" "b")

    for inst in ${inst_list[@]};
    do
        source_num=`ls ${inst}/x*.reg | wc -l`
        source_index=1
        while [ $source_index -le ${source_num} ]
        do
            mkdir x${source_index}
            cd x${source_index}

            evt_file=`ls ../fxt/products/*fxt_${inst}*.fits`
            echo "

read events ${evt_file}
.
yes
filter region ../${inst}/x${source_index}.reg
extract spectrum
save spectrum x${source_index}_${inst}.pi clobber=yes
clear region
filter region ../bkg.reg
extract spectrum
save spectrum bkg_${inst}.pi clobber=yes
exit
no" > xselect_para_${inst}.xco

            xselect @xselect_para_${inst}.xco
            fxtarfgen specfile=x${source_index}_${inst}.pi expfile=../${inst}/expo_${inst}-without_vign.img outfile=x${source_index}_${inst}.arf clobber=yes
            fxtrmfgen specfile=x${source_index}_${inst}.pi outfile=x${source_index}_${inst}.rmf clobber=yes
            cd ..
            let source_index++
        done
    done
    cd ../..
done
