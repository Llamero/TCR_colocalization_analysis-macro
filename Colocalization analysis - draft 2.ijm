setBatchMode(true);
run("Close All");
//run("Clear Results");
use_mask = true;
sample_dir = getDirectory("Choose input directory");
file_list = getFileList(sample_dir);
sample_list = newArray(0);
bead_list = newArray(0);
error_string = "";
for(a=0; a<file_list.length; a++){
	if(matches(file_list[a], "^(?!.+Bead|.+bead).+\.czi$")){ //Search for all czi files that do not contain the word "bead"
		test_string = replace(file_list[a], "\\)-", "\\)-Bead-");
		test_string = replace(file_list[a], "\\)_", "\\)-Bead_");
		test_string = replace(test_string, "-[0-9]+\\.czi$", "");
		for(b=0; b<file_list.length; b++){
			if(startsWith(file_list[b], test_string)){ //Search for tbe bead control image that corresponds to the sample image
				sample_list = Array.concat(sample_list, file_list[a]);
				bead_list = Array.concat(bead_list, file_list[b]);
				break;
			}
		}
		if(b == file_list.length) error_string += "Bead image not found for following stack: " + file_list[a] + "\r\n";
	}
}

for(a=0; a<sample_list.length; a++){
	run("Bio-Formats Importer", "open=[" + sample_dir + sample_list[a] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	sample = getTitle();
	run("Bio-Formats Importer", "open=[" + sample_dir + bead_list[a] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	bead = getTitle();
	processImages(bead, sample, a, sample_list[a]);
	measureColocalization("C1-registered time points", "C2-registered time points", a);
	run("Close All");
}
saveAs("Results", sample_dir + "Results.csv");
selectWindow("Results");
print("\\Clear");
print(error_string);
setBatchMode("exit and display");

function measureColocalization(bead, sample, number){
	for(c=0; c<2; c++){
		if(c == 0){
			id = "Bead";
			stack = bead; 
		}
		else{
			id = "Sample";
			stack = sample;
		}
		selectWindow(stack);
		getDimensions(width, height, i_channels, i_slices, i_frames);
		run("Split Channels");
		merge_string = "";
		if(use_mask){
			for(a=1; a<=i_channels; a++){
					selectWindow("C" + a + "-" + stack);
					run("Duplicate...", "title=[Mask " + a + "] duplicate");
					selectWindow("Mask " + a);
					setAutoThreshold("Triangle dark stack");
					run("Convert to Mask", "method=Triangle background=Dark");
					merge_string += "c" + a + "=[Mask " + a + "] ";
			}
			run("Merge Channels...", merge_string + "create keep ignore");
			selectWindow("Composite");
			saveAs("tiff", sample_dir + "Mask " + number);
			close("Mask " + number + ".tif");
		}
	
		for(a=1; a<i_channels; a++){
			for(b=a+1; b<=i_channels; b++){
				print("\\Clear");
				wait(200);
				if(use_mask){
					imageCalculator("Add create stack", "Mask " + a,"Mask " + b);					
					run("Coloc 2", "channel_1=[C" + a + "-" + stack + "] channel_2=[C" + b + "-" + stack + "] roi_or_mask=[Result of Mask " + a + "] threshold_regression=Bisection li_histogram_channel_1 li_histogram_channel_2 li_icq spearman's_rank_correlation manders'_correlation kendall's_tau_rank_correlation 2d_intensity_histogram costes'_significance_test psf=4 costes_randomisations=10");
				}
				else run("Coloc 2", "channel_1=[C" + a + "-" + stack + "] channel_2=[C" + b + "-" + stack + "] roi_or_mask=<None> threshold_regression=Bisection li_histogram_channel_1 li_histogram_channel_2 li_icq spearman's_rank_correlation manders'_correlation kendall's_tau_rank_correlation 2d_intensity_histogram costes'_significance_test psf=4 costes_randomisations=10");
				close("Result of Mask " + a);
				wait(200);
				log_string = getInfo("log");
				pearson_end_index = indexOf(log_string, "Pearson's R value (below threshold),");
				pearson_t_end_index = indexOf(log_string, "Li's ICQ value,");
				spearman_end_index = indexOf(log_string, "Spearman's correlation t-statistic, ");
				pearson_coefficient = parseFloat(substring(log_string, pearson_end_index-5,  pearson_end_index-1));
				pearson_t_coefficient = parseFloat(substring(log_string, pearson_t_end_index-5,  pearson_t_end_index-1));
				spearman_coefficient = parseFloat(substring(log_string, spearman_end_index-11,  spearman_end_index-1));
				setResult(id + " C" + a + "-C" + b + " Pearson", nResults-1, pearson_coefficient);
				setResult(id + " C" + a + "-C" + b + " Pearson (above threshold)", nResults-1, pearson_t_coefficient);
				setResult(id + " C" + a + "-C" + b + " Spearman", nResults-1, spearman_coefficient);
				updateResults();
			}	
		}
		for(a=1; a<=i_channels; a++) close("Mask " + a);
	}
}


function processImages(i, j, number, title){
	selectWindow(i);
	getDimensions(width, height, i_channels, i_slices, i_frames);
	selectWindow(j);
	getDimensions(dummy, dummy, j_channels, j_slices, j_frames);
	
	//Make stacks equal # of slices
	selectWindow(i);
	while(i_slices < j_slices){
		Stack.setSlice(i_slices);
		run("Add Slice", "add=slice");
		getDimensions(width, height, i_channels, i_slices, i_frames);
	}
	selectWindow(j);
	while(j_slices < i_slices){
		Stack.setSlice(j_slices);
		run("Add Slice", "add=slice");
		getDimensions(width, height, j_channels, j_slices, j_frames);
	}
	
	//Remove channel 4 from analysis
	selectWindow(i);
	run("Split Channels");
	close("C4-" + i);
	run("Merge Channels...", "c1=[C1-" + i + "] c2=[C2-" + i + "] c3=[C3-" + i + "] create");
	selectWindow(j);
	run("Split Channels");
	close("C4-" + j);
	run("Merge Channels...", "c1=[C1-" + j + "] c2=[C2-" + j + "] c3=[C3-" + j + "] create");
	
	//remove border
	x_min = -1;
	x_max = -1;
	y_min = -1;
	y_max = -1;
	selectWindow(i);
	Stack.setChannel(1);
	for(a=0; a<width; a++){
		makeRectangle(a, 0, 1, height);
		getStatistics(dummy, dummy, dummy, max);
		if(max > 0 && x_min < 0) x_min = a;
		if(max == 0 && x_min >= 0){
			x_max = a;
			break;
		}
	}
	for(a=0; a<height; a++){
		makeRectangle(0, a, width, 1);
		getStatistics(dummy, dummy, dummy, max);
		if(max > 0 && y_min < 0) y_min = a;
		if(max == 0 && y_min >= 0){
			y_max = a;
			break;
		}
	}
	makeRectangle(x_min, y_min, x_max-x_min, y_max-y_min);
	run("Crop");
	selectWindow(j);
	makeRectangle(x_min, y_min, x_max-x_min, y_max-y_min);
	run("Crop");
	
	//register channels - split - remove padding - save copy of stack
	print("\\Clear");
	selectWindow(i);
	run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
	selectWindow(j);
	run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
	run("Merge Channels...", "c1=[" + i + "] c2=[" + j + "] create");
	wait(200);
	selectWindow("Merged");
	run("Correct 3D drift", "channel=1 multi_time_scale sub_pixel only=0 lowest=1 highest=" + i_slices);
	selectWindow("registered time points");
	run("Split Channels");
	selectWindow("C1-registered time points");
	run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
	slice_max = 0;
	while(slice_max == 0){
		Stack.setSlice(nSlices);
		getStatistics(dummy, dummy, dummy, slice_max);
		if(slice_max == 0) run("Delete Slice", "delete=slice");
	}
	saveAs("tiff", sample_dir + "Registered Beads " + number);
	rename("C1-registered time points");
	selectWindow("C2-registered time points");
	run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
	slice_max = 0;
	while(slice_max == 0){
		Stack.setSlice(nSlices);
		getStatistics(dummy, dummy, dummy, slice_max);
		if(slice_max == 0) run("Delete Slice", "delete=slice");
	}
	saveAs("tiff", sample_dir + "Registered Sample " + number);
	rename("C2-registered time points");
	close("Merged");
	setBatchMode("exit and display");
	exit();
	//Print shift vectors to results table
	setResult("Sample ID", nResults, title);
	shift_array = newArray(i_channels-1);
	log_string = getInfo("log");
	for(a=i_channels-2; a>=0; a--){
		start_index = indexOf(log_string, "    frame " + a+1 + " correcting drift ");
		shift_array[a] = substring(log_string, start_index+28, lengthOf(log_string)-1);
		log_string = substring(log_string, 0, start_index);
	}
	Array.print(shift_array);
	for(a=0; a<shift_array.length; a++){
		setResult("Channel " + a+1 + " shift", nResults-1, shift_array[a]);
	}
	updateResults();
}

