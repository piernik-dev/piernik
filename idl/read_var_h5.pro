FUNCTION READ_VAR_H5, file=file, var=var, start=start, count=count

   file_id = H5F_OPEN(file)

      dst_id = H5D_OPEN(file_id, var)
      dsp_id = H5D_GET_SPACE(dst_id)

      H5S_SELECT_HYPERSLAB, dsp_id, start, count,/RESET
      mem_id = H5S_CREATE_SIMPLE(count)

      data_arr  = H5D_READ(dst_id,FILE_SPACE=dsp_id, MEMORY_SPACE=mem_id)

   H5S_CLOSE, mem_id
   H5S_CLOSE, dsp_id
   H5D_CLOSE, dst_id
   H5F_CLOSE, file_id

   RETURN, data_arr

END
