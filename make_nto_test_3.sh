#!/bin/bash
# Script to run NTOs
#NTOs calculation script from the TDDFT chk file and a list state file name nto_states 
#Scrip Revised by Mohammed A Jabed, NDSU 
##script revised by Sean Ferrell sferr092@fiu.edu, improved speed so as to run ~2-4x faster, removed redundancy.
##to estimate time to run a decent approximate is 7.6*10^-5 *(n^2) = hours, where n is the number of atoms.


file=$(rev <<< "$1" | cut -d"." -f2- | rev)
file_nto="${file}_nto"

mkdir -p ntos
mkdir -p ntos/tmp/
rm ntos/tmp/* 

module load gaussian/16-c02-avx2
export GAUSS_MDEF=20GB
export GAUSS_MEMDEF=20GB
export GAUSS_PDEF=32
export GAUSS_UDEF=20GB

echo "Copying $1 and NTO templates to script"
cp "$1" "ntos/$file_nto.chk" ###if rerunning comment this out
cp nto_states ntos/nto_states.nto

cd ntos || exit 

#formchk "$file_nto".chk "${file}".fchk ###if rerunning comment this line out
 
a1=$(grep -n Route  "${file}".fchk  | head -1 | awk -F':' '{print $1}') 
a=$(echo "$a1" +1 | bc -l ) 
b1=$(grep -n Charge  "${file}".fchk  | head -1 | awk -F':' '{print $1}') 
b=$(echo "$b1" -1 | bc -l ) 

Route=$(head -"$b" "$file".fchk | tail -n +"$a"  | paste -sd "" )  
echo $Route
solvent=$(echo "$Route" |  grep -oi 'scrf\S*' )
echo $solvent
pseudo=$(echo "$Route" |  grep -oi 'pseudo\S*' ) 
functional=$(head -2 "${file}".fchk | tail -1  | awk '{print $2 }')


echo "%mem=16GB
%nprocshared=32
%chk=file_nto1.chk

#p guess=(read,only) density=(check,transition=nto_trans) nosymm $pseudo $solvent geom=Allcheck 
#p  pop=(NTO,saveNTO) ${functional}/chkbasis  
" > "${file}"_nto1.com  


echo "%mem=16GB
%nprocshared=32
%chk=file_nto2.chk

#p guess=(read,only) density=(check,transition=nto_trans) nosymm $pseudo  $solvent geom=AllCheck 
#p gfprint pop=(minimal,NTO,SaveNTO) ${functional}/chkbasis 

" > "${file}"_nto2.com  


echo "Generating .com files and running Gaussian"
#source nto_states.nto

for state in $(cat nto_states.nto ) 
do
 echo "Running Gaussian for ${file_nto}1_state${state}.com"
 cp "${file}_nto1.com" "${file_nto}1_state${state}.com"
 cp "${file}_nto2.com" "${file_nto}2_state${state}.com"
 cp "${file}_nto.chk" "${file_nto}1_state${state}.chk" 
 sed -i "s/nto_trans/$state/g" "${file_nto}1_state${state}.com"
 sed -i "s/nto_trans/$state/g" "${file_nto}2_state${state}.com"
 sed -i "s/file_nto1/${file_nto}1_state${state}/g" "${file_nto}1_state${state}.com"
 sed -i "s/file_nto2/${file_nto}2_state${state}/g" "${file_nto}2_state${state}.com"
 
 export G16_INPUT="${file_nto}1_state${state}.com"
 export G16_OUTPUT="${file_nto}1_state${state}.log"
 #export GAUSS_SCRDIR=/mmfs1/projects/svetlana.kilina/g16.$PBS_JOBID
 export GAUSS_SCRDIR=$SCRATCH/g16.$PBS_JOBID
 mkdir -p "$GAUSS_SCRDIR"

 g16 < "$G16_INPUT" > "$G16_OUTPUT"

 if grep -q "Normal term" "$G16_OUTPUT"; then 
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
 	echo " NTO1 calculation of the state $state is finished normally" 
 	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
 else 
 	echo "!!!!!!!!"
 	echo "failure to do NTO1 $state"
	echo "!!!!!!"
 fi 
 
 rm -rf "$GAUSS_SCRDIR"
 echo "Creating  ${file_nto}2_state${state}.chk"
 cp "${file_nto}1_state${state}.chk" "${file_nto}2_state${state}.chk" 
 
 echo "Running Gaussian for ${file_nto}2_state${state}.com"
 export G16_INPUT="${file_nto}2_state${state}.com"
 export G16_OUTPUT="${file_nto}2_state${state}.log"
 #export GAUSS_SCRDIR=/gpfs1/projects/svetlana.kilina/g16.$PBS_JOBID
 export GAUSS_SCRDIR=$SCRATCH/g16.$PBS_JOBID
 mkdir -p "$GAUSS_SCRDIR"

 g16 < "$G16_INPUT" > "$G16_OUTPUT"

 if grep -q "Normal term" "$G16_OUTPUT"; then 
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	echo " NTO2 calculation of the state $state finish normally" 
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
 else 
 	echo "!!!!!!!!"
 	echo "failure to do NTO2 $state"
	echo "!!!!!!"

 fi 
  
 rm -rf "$GAUSS_SCRDIR"
  
 formchk -3 "${file_nto}2_state${state}.chk"
 echo "got to NTO2_$state formchk" $(date) >> status
##################################################
alpha=$(grep 'alpha electrons' "${file}"_nto2_state"${state}".log | awk '{print $1}' )

grep 'Alpha  occ. ' "${file}"_nto2_state"${state}".log  | cut -c32-85 > tmp/OCC1_"${file}"_"${state}".log 
grep -m2 'Alpha virt' "${file}"_nto1_state"${state}".log  | cut -c32-85 > tmp/virt1_"${file}"_"${state}".log
for i in $(cat tmp/OCC1_"${file}"_"${state}".log) ;  do echo "$i" >> tmp/OCC2_"${file}"_"${state}".log ; done
for i in $(cat tmp/virt1_"${file}"_"${state}".log) ;  do echo "$i" >> tmp/virt2_"${file}"_"${state}".log ; done

a=$(nl -b a tmp/OCC2_"${file}"_"${state}".log \
  | awk '$2>=.2' \
  | cut -c3-6 \
  | head -1 \
  | sed 's/$/+1/' | bc -l) ##previously did the occupied and then unoccupied which should be wrong
#$2>=.2 is twenty percent, change .2 in both spots to change percentage

a_virt=$(echo "$a-1" | bc -l) 
#$2>=.2 is twenty percent, change .2 in both spots to change percentage

n=$(nl -b a tmp/OCC2_"${file}"_"${state}".log \
  | awk '$2>=.2' \
  | wc \
  | awk '{print $1}') 
#$2>=.2 is twenty percent, change .2 in both spots to change percentage

n_virt=$(nl -b a tmp/virt2_"${file}"_"${state}".log \
  | awk '$2>=.2' \
  | wc \
  | awk '{print $1}') 
#$2>=.2 is twenty percent, change .2 in both spots to change percentage

echo "$n" 
b=$(echo "$alpha+$n" | bc -l )
echo "$a" 
echo "$b" 
echo $(seq $a 1 $b) "hole" >> status

echo "$n_virt" 
b_virt=$(echo "$alpha-$n_virt+1" | bc -l )
echo "$(echo "$a_virt")" 
echo "$b_virt" 
echo $(seq $b_virt 1 $a_virt) "electron" >> status

##setting up functions to split the cubegen task to increase speed. max effective cpus for cubegen was found to be 16.
##these functions can be called in parallel to enable a 2-4x speed up. further splitting is possible but number of occupied holes or electrons is normally =<2
holeodds() {
 for i in $(seq "$a" 2 "$b") ; do
  echo "cubegen mo="$i"" $(date) >> status
  cubegen 16 mo="$i" "${file}_nto2_state${state}.fchk" "${file}_nto1_state${state}_$i.cube" -2 h
 done 
}
holeevens() {
 for ii in $(seq $(($a+1)) 2 "$b") ; do
  echo "cubegen mo="$ii"" $(date) >> status
  cubegen 16 mo="$ii" "${file}_nto2_state${state}.fchk" "${file}_nto1_state${state}_$ii.cube" -2 h
 done
}
elecevens() {
 for x in $(seq "$b_virt" 2 "$a_virt") ; do
  echo "cubegen mo="$x"" $(date) >> status
  cubegen 16 mo="$x" "${file}_nto2_state${state}.fchk" "${file}_nto2_state${state}_$x.cube" -2 h
 done
}
elecodds() {
 for xx in $(seq $(($b_virt+1)) 2 "$a_virt") ; do
  echo "cubegen mo="$xx"" $(date) >> status
  cubegen 16 mo="$xx" "${file}_nto2_state${state}.fchk" "${file}_nto2_state${state}_$xx.cube" -2 h
 done
}
##calling functions simultaneously 
holeodds &
holeevens &
elecodds &
elecevens & 

################################################
done
rm -r tmp/
wait
echo "succesfully generated NTO cube files"
j=$(LC_ALL=C date +%I+%M+%S)
printf "%s\n" "$j"

