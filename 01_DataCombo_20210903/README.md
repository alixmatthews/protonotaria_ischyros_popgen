# Combining data  
Alix Matthews

Date: 20210903

## Computer and path of working directory
- AHPCC: ```/scrfs/storage/amatthews/20210816_projects/20210816_snp/01_DataCombo_20210903```

## Description of what you're doing and why
- Combining previous data (from 20210412_snp) with new data so they are all in the same directory

## Code used to do so

```
cd /scrfs/storage/amatthews/20210816_projects/20210816_snp

mkdir 01_DataCombo_20210903

cd 01_DataCombo_20210903

mkdir Adapter_Removed_bb
``` 

```
cp -R /scrfs/storage/amatthews/20210412_snp/00_PP_20210427/Adapter_Removed_bb/* /scrfs/storage/amatthews/20210816_projects/20210816_snp/01_DataCombo_20210903/Adapter_Removed_bb
```



```
cp -R /scrfs/storage/amatthews/20210816_projects/20210816_snp/00_PP_20210902/Adapter_Removed_bb/* /scrfs/storage/amatthews/20210816_projects/20210816_snp/01_DataCombo_20210903/Adapter_Removed_bb
```







