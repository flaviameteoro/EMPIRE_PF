program gcom_test
use mpl
implicit none
integer  nproc_um_npes,nproc_max, e_um_npes,err,mype_id, mycolour,mype
integer int_log,n_procs,info,couple_colour,couple_myrank,couple_nproc
integer couple_root


      INTEGER MY_MPI_COMMUNICATOR      ! Local MPI Communicator
      INTEGER GLOBAL_MPI_COMMUNICATOR  ! Global MPI Communicator
      integer COUPLE_MPI_COMMUNICATOR

      LOGICAL IS_UMCET_ENSEMBLE        ! Indicates if ensemble
      integer MPL_UNDEFINED
      MPL_UNDEFINED=mpi_undefined
      nproc_max=1024
      nproc_um_npes = nproc_max
!      CALL GC_INIT(dummy_env,mype,nproc_max)
! added lines for UMCET
! Check to see if we are part of an UMCET ensemble

       e_um_npes=32
        CALL GC_INIT_INTRO(GLOBAL_MPI_COMMUNICATOR)
        CALL MPL_COMM_RANK(GLOBAL_MPI_COMMUNICATOR,mype_id,err)
!       Using division of integers to get the "colours" for the
!       split. This make sure that the ensemble members are
!       close to each other on their origional MPI ID
!       and probably close together in the machine too
        mycolour = mype_id/ e_um_npes
!        write(6,*)'hello1',mype_id,mycolour
        CALL MPL_COMM_SPLIT(GLOBAL_MPI_COMMUNICATOR,mycolour, &
                           mype_id,MY_MPI_COMMUNICATOR,err)

        CALL GC_INIT_FINAL(mype,nproc_max,MY_MPI_COMMUNICATOR)

        if (mype.eq.0) then
           couple_colour=9999
        else
           couple_colour=MPL_UNDEFINED
        endif

        CALL MPL_COMM_SPLIT(GLOBAL_MPI_COMMUNICATOR,couple_colour, &
                           mype_id,COUPLE_MPI_COMMUNICATOR,err)

        if (mype.eq.0) then 
     call MPl_Comm_Rank (COUPLE_MPI_COMMUNICATOR, couple_myRank, err)
     call MPl_Comm_Size (COUPLE_MPI_COMMUNICATOR, couple_nProc, err)

         call MPL_allreduce(0,couple_root,1,MPL_INTEGER8,MPL_MAX,&
             COUPLE_MPI_COMMUNICATOR, err)
        write(6,"(a,i2,i2,i2,i2,i2,i2)") 'hello2',mype_id,mycolour,mype,couple_myRank,couple_nProc,couple_root
        endif
int_log=mype_id
CALL GC_IMAX(1,e_um_npes,info,int_log)
        write(6,"(a,i2,i2,i2,i2)") 'hello3',mype_id,mycolour,mype,int_log
        CALL GC_EXIT()

        end program

