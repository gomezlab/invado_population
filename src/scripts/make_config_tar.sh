cd ../../data; 
tar cf all_configs.tar * --exclude *.mp4 --exclude *.mov --exclude *.avi* --exclude *.tiff --exclude *.TIF --exclude *.tif --exclude *.png --exclude config;
if [ -n "$1" ]; then
	scp all_configs.tar $1:
	rm all_configs.tar
fi
