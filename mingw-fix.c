#if defined(__MINGW32__) || defined(__MINGW64__)
/*double inline drand48(int seed) {
	return 1.0*rand();
}
double inline srand48(int seed) {
	return 1.0*rand();
}*/
double inline drand48() {
	return 1.0*rand();
}
double inline srand48() {
	return 1.0*rand();
}
#endif