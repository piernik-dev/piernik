find . -name "*.F90" -type f -exec \
sed -e "s/diff_nml(\(.*\))/if (.not.nh%initialized) call nh%init()\n\
         open(newunit=nh%lun, file=nh%tmp1, status=\"unknown\")\n\
         write(nh%lun,nml=\1)\n\
         close(nh%lun)\n\
         open(newunit=nh%lun, file=nh%par_file)\n\
         nh%errstr=\"\"\n\
         read(unit=nh%lun, nml=\1, iostat=nh%ierrh, iomsg=nh%errstr)\n\
         close(nh%lun)\n\
         call nh%namelist_errh(nh%ierrh, \"\1\")\n\
         read(nh%cmdl_nml,nml=\1, iostat=nh%ierrh)\n\
         call nh%namelist_errh(nh%ierrh, \"\1\", .true.)\n\
         open(newunit=nh%lun, file=nh%tmp2, status=\"unknown\")\n\
         write(nh%lun,nml=\1)\n\
         close(nh%lun)\n\
         call nh%compare_namelist()/g" -i {} \;
