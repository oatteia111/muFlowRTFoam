	int mpi_tasks,mpi_myself;
	#ifdef USE_MPI
		MP_TYPE comm = MPI_COMM_WORLD;
		if (MPI_Init(&argc, &argv) != MPI_SUCCESS) { return EXIT_FAILURE; }
		if (MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks) != MPI_SUCCESS) {return EXIT_FAILURE;}
		if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS){exit(4);}
		std::cerr << "MPI started. " << mpi_myself << std::endl;
		if (mpi_myself == 0)
		{
			freak.setDB (cur_dir/"phreeqc.dat");
			freak.setChemFile(cur_dir/"initChem.pqi");
			freak.setData(ph_data);
			nxyz=ph_data[0];
		}
		MPI_Bcast(&nxyz, 1, MPI_INT, 0, MPI_COMM_WORLD);
		PhreeqcRM phreeqc_rm(nxyz, MPI_COMM_WORLD);
		freak.PhreeqcRM_ptr = &phreeqc_rm;
		if (mpi_myself > 0)
		{
			std::cerr << "Started MPI worker " << mpi_myself << std::endl;
			phreeqc_rm.MpiWorker();
			//MPI_Finalize(); // does not seem to change anything
			return EXIT_SUCCESS;
		}
	#endif
