import java.util.Arrays;

public class StepCounter {
	private static final int WINDOW_LENGTH = 5; //"range" for noise smoothing
	private static final double DEVIATION_SCALAR= 0.8; //constant that is multiplied with (mean+standard deviation) to get threshold value for naive algorithm
	private static final int THRESHOLD_MULTIPLE=5; //multiplication factor for threshold to ensure it is above the magnitudes of invalid steps
	private static final double MEAN_SCALAR= 0.54; //scalar applied to mean value of data to calculate minimum reasonable data value (anything smaller in mag. cannot be a step)
	private static final double MEAN_SHIFT=0.72; //value to shift mean value of data from to calculate minimum reasonable data value (anything smaller in mag. cannot be a step)
	private static final double PERCENT_ACCEPTANCE_POTENTIAL_STEP_ACTIVITY=0.8;
	public static void main(String[] args) {
		String[] columnNames={"time", "accel-x", "accel-y", "accel-z"};
		CSVData test = new CSVData("data/50StepWalkFemale(2).csv", columnNames, 1);
		test.correctTime(test);
		System.out.println("Improved Algorithm:"+countSteps(test.getColumn(0),test.getRows(1,test.getNumRows()-1), 10));
		System.out.println("Old Algorithm:"+countSteps(test.getColumn(0),test.getRows(1,test.getNumRows()-1)));
	}
	/***
	 * NOTE: IMPROVED ALGORITHM 
	 * Counts the number of steps with an adaptive threshold that is calculated for each sensorData value. 
	 * If the resultant of all 3 components of acceleration (from sensorData) is above the adaptive threshold for that point
	 * and greater in magnitude than its two adjacent points, it is considered a step.
	 * @param times array of times corresponding to each resultant magnitude from the 3 components of acceleration from sensorData
	 * @param sensorData 2-D array of sensor data including times(ms) and the three acceleration components (x,y,z)
	 * @param windowLength length of window adjacent to each resultant acceleration that the thresholds are calculated from
	 * 		Ex: windowLength=3 at times=5 means thresholds will be calculated from resultant accelerations corresponding to times
	 * 		[2,8]
	 * @return number of steps counted from sensorData
	 */
	private static int countSteps(double[] times, double[][] sensorData, int windowLength) {
		int stepCount = 0;
		double[] arr = new double[times.length];
		arr = calculateMagnitudes(sensorData);
		double mean = calculateMean(arr);
		double minLimit = calculateMinLimit(mean); //minimum magnitude to be a reasonable step
//		arr = NoiseSmoothing.generalRunningAverage(arr, 10);
		double[] thresholds = calculateAdaptiveThreshold(arr, windowLength, minLimit);
		for(int i = 1; i < arr.length-1; i++) {
			if (arr[i] > arr[i-1] && arr[i] > arr[i+1]) {
				if(arr[i] > thresholds[i]) {
					stepCount++;
//					System.out.println(stepCount +" " + times[i]);
				}
			}
		}
		return stepCount;	
	}	
	/**
	 * NOTE: OLD NAIVE ALGORITHM
	 * Counts the number of steps with a static threshold that is calculated by adding the mean of sensorData resultant acceleration
	 * vector with the standard deviation of sensorData resultant acceleration vector and multiplied by DEVIATION_SCALAR
	 * @param times array of times corresponding to each resultant magnitude from the 3 components of acceleration from sensorData
	 * @param sensorData 2-D array of sensor data including times(ms) and the three acceleration components (x,y,z)
	 * @return number of steps counted from sensorData
	 */
	private static int countSteps(double[] times, double[][] sensorData) {
		int stepCount = 0;
		double[] arr = new double[times.length];
		arr = calculateMagnitudes(sensorData);
		double mean = calculateMean(arr);
		double deviation = calculateStandardDeviation(arr, mean);
//		System.out.println(mean+deviation);
		for(int i = 1; i < arr.length-1; i++) {
			if (arr[i] > arr[i-1] && arr[i] > arr[i+1]) {
				if(arr[i] > (mean+deviation)*DEVIATION_SCALAR) {
					stepCount++;
//					System.out.println(stepCount +" " + times[i]/1000);
				}
			}
		}
		
		return stepCount;
	}
	/**
	 * Returns minimum magnitude of acceleration resultant vector to be a reasonable step. 
	 * MEAN_SCALAR and MEAN_SHIFT derived from trial-and-error of testing for the best minLimit for different step situations,
	 * and finding the best linear fit for the data
	 * @param mean the average from the array of acceleration resultant vectors
	 * @return minimum magnitude of acceleraiton resultant vecotr to be a reasonable step
	 */
	public static double calculateMinLimit(double mean){
		return MEAN_SCALAR*mean+MEAN_SHIFT;
	}
	/**
	 * Returns array of threshold values for input window size away from each arr value. If magnitude is unreasonably 
	 * low(not a step), threshold will be set to an extreme high value to avoid counting peaks in no-step noise areas.
	 * @param arr array of magnitudes to calculate thresholds from
	 * @param windowLength length of "range" of values next to each magnitude value to calculate threshold from
	 * @return array of threshold values for every magnitude value calculated with int windowLength values to the left/right of
	 * magnitude value
	 */
	public static double[] calculateAdaptiveThreshold(double[] arr, int windowLength, double minLimit) {
		double[] result=new double[arr.length];
		for (int i = 0; i < arr.length; i++){
			if(i+windowLength < arr.length) { //check if in right-hand range
				if(i-windowLength >= 0) { //check if in left-hand range
					double meanForInterval = calculateMeanInInterval(arr, i-windowLength, i+windowLength);
					double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, i-windowLength, i+windowLength);
					if(isMagnitudeUnreasonablySmall(arr, WINDOW_LENGTH, i, minLimit))
						result[i]=((meanForInterval+deviationForInterval))*THRESHOLD_MULTIPLE;
					else
						result[i] = (meanForInterval+deviationForInterval);
				} else { //if index is within windowLength from the left-hand end of array
					double meanForInterval = calculateMeanInInterval(arr, 0, i+windowLength);
					double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, 0, i+windowLength);
					if(isMagnitudeUnreasonablySmall(arr, WINDOW_LENGTH, i+WINDOW_LENGTH, minLimit)){
						result[i]=((meanForInterval+deviationForInterval))*THRESHOLD_MULTIPLE;
					}else{
						result[i] = (meanForInterval+deviationForInterval);
					}
				}
			} else {//if index is within windowLength from the right-hand end of array
				double meanForInterval = calculateMeanInInterval(arr, i-windowLength, arr.length);
				double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, i-windowLength, arr.length);
				if(isMagnitudeUnreasonablySmall(arr, WINDOW_LENGTH, i-WINDOW_LENGTH, minLimit)){
					result[i]=((meanForInterval+deviationForInterval))*THRESHOLD_MULTIPLE;
				}else{
					result[i] = (meanForInterval+deviationForInterval);
				}
			}
		}
		return result;
	}
	
//	public static double[] calculateWindow(double[] arr, int windowLength) {
//		double[] result=new double[arr.length];
//		for (int i = 0; i < arr.length; i++){
//			if(i+windowLength < arr.length) {
//				if(i-windowLength >= 0) {
//					double meanForInterval = calculateMeanInInterval(arr, i-windowLength, i+windowLength);
//					double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, i-windowLength, i+windowLength);
//					result[i] = (meanForInterval+deviationForInterval);
//				} else {
//					double meanForInterval = calculateMeanInInterval(arr, 0, i+windowLength);
//					double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, 0, i+windowLength);
//					result[i] = (meanForInterval+deviationForInterval);
//				}
//			} else {
//				double meanForInterval = calculateMeanInInterval(arr, i-windowLength, arr.length);
//				double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, i-windowLength, arr.length);
//				result[i] = (meanForInterval+deviationForInterval);
//			}
//		}
//		return result;
//	}
	
	/**
	 * Checks if the magnitudes of arr are unreasonably small (insignificant data that has magnitudes too small to count
	 * as steps) within a windowLength range from the arr value index
	 * @param arr magnitude(accel) values to check from
	 * @param windowLength range from the arr value
	 * @param index arr index value to check magnitudes from
	 * @param threshold threshold for reasonable magnitudes (anything under is considered too small to count as possible steps)
	 * @return true if >80% of values are under the threshold (aka too small to be considered any significant step activity)
	 * false if <80% of values are under the threshold (aka of enough significance to be considered for potential step activity)
	 */
	private static boolean isMagnitudeUnreasonablySmall(double[] arr, int windowLength, int index, double threshold) {
		int underThresholdCounter = 0, totalIterations=0;
		for(int i = index-windowLength; i < index+windowLength; i++) {
			if(arr[i] < threshold) 
				underThresholdCounter++;
			totalIterations++;
		}
		if((double)underThresholdCounter/(double)totalIterations > PERCENT_ACCEPTANCE_POTENTIAL_STEP_ACTIVITY) 
			return true;
		return false;
	}
//	/**
//	 * 
//	 * @param magnitudes
//	 * @param averageLength
//	 * @return
//	 */
//	public static double[] noiseSmoothing(double[] magnitudes, int averageLength) {
//		double[] result = new double[magnitudes.length-averageLength];
//		result = NoiseSmoothing.generalRunningAverage(magnitudes, averageLength);
//		return result;
//	}
	/**
	 * Returns the magnitude of the resultant vector of three input magnitudes of vectors
	 * @param x magnitude of x-component vector
	 * @param y magnitude of y-component vector
	 * @param z magnitude of z-component vector
	 * @return magnitude of the resultant vector from three input component vectors
	 */
	public static double calculateMagnitude(double x, double y, double z) {
		return Math.sqrt(x*x + y*y + z*z);
	}
	/**
	 * Returns a 1D array with magnitudes of the resultant vectors of three vector components for multiple rows
	 * Calculates the magnitude of each row of three vector components for multiple rows of sensorData
	 * @param sensorData 2-D array of sensor data including times(ms) and the three acceleration components (x,y,z)
	 * @return 1D array with magnitudes of the resultant vectors of three vector components for multiple rows
	 */
	private static double[] calculateMagnitudes(double[][] sensorData) {
		double[] result = new double[sensorData.length];
		for(int i = 0; i < sensorData.length; i++) {
			double x = sensorData[i][1];
			double y = sensorData[i][2];
			double z = sensorData[i][3];
			result[i] = calculateMagnitude(x,y,z);
		}
		return result;
	}
	/**
	 * Returns standard deviation of a 1D array of values
	 * @param arr 1D array of values to calculate standard deviation from
	 * @param mean average of the 1D array 
	 * @return standard deviation of arr
	 */
	private static double calculateStandardDeviation(double[] arr, double mean) {
		double sum = 0;
		for (int i = 0; i < arr.length; i++) {
			sum += (arr[i] - mean)*(arr[i] - mean);
		}
		return Math.sqrt(sum/(arr.length-1));
	}
	/**
	 * Returns standard deviation of a 1D array in a select interval
	 * @param arr 1D array of values to calculate standard deviation from
	 * @param mean average of the 1D array 
	 * @param startInterval starting index of interval
	 * @param endInterval ending index of interval (non-inclusive)
	 * @return standard deviation of arr in interval
	 */
	private static double calculateStandardDeviationInInterval(double[] arr, double mean, int startInterval, int endInterval) {
		double sum = 0;
		for (int i = startInterval; i < endInterval; i++) {
			sum += (arr[i] - mean)*(arr[i] - mean);
		}
		return Math.sqrt(sum/(double)(endInterval-startInterval));
	}
	/**
	 * Returns average from 1D array of values
	 * @param arr 1D array of values
	 * @return average of 1D array 
	 */
	private static double calculateMean(double[] arr) {
		double sum = 0;
		for(int i = 0; i < arr.length; i++) {
			sum += arr[i];
		}
		return sum/arr.length;
	}
	/**
	 * Returns average of 1D array in a select interval
	 * @param arr 1D array of values
	 * @param startIndex starting index of interval
	 * @param endIndex ending index of interval (non-inclusive)
	 * @return average of 1D array
	 */
	private static double calculateMeanInInterval(double[] arr, int startIndex, int endIndex) {
		double sum = 0;
		for(int i = startIndex; i < endIndex; i++) {
			sum += arr[i];
		}
		return sum/(double)(endIndex-startIndex);
	}
}
