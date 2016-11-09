import java.util.Arrays;

public class StepCounter {
	private static final int WINDOW_LENGTH = 5;
	private static final double DEVIATION_SCALAR= 0.8;
	private static final int THRESHOLD_MULTIPLE=5; //multiplication factor for threshold to ensure it is above the magnitudes of invalid steps
	public static void main(String[] args) {
		String[] columnNames={"time", "gyro-x", "gyro-y", "gyro-z"};
		CSVData test = new CSVData("data/50StepWalkFemale(2).csv", columnNames, 1);
		test.correctTime(test);
		System.out.println("Improved Algorithm:"+countSteps(test.getColumn(0),test.getRows(1,test.getNumRows()-1), 10));
		System.out.println("Old Algorithm:"+countSteps(test.getColumn(0),test.getRows(1,test.getNumRows()-1)));
	}
	//improved algorithm
	private static int countSteps(double[] times, double[][] sensorData, int windowLength) {
		int stepCount = 0;
		double[] arr = new double[times.length];
		arr = calculateMagnitudes(sensorData);
		double mean = calculateMean(arr);
		double minLimit = calculateMinLimit(mean); //minimum magnitude to be a reasonable step
//		arr = NoiseSmoothing.generalRunningAverage(arr, 10);
		double[] thresholds = calculateWindow(arr, windowLength, minLimit);
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
	//old algorithm
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
	
	public static double calculateMinLimit(double mean){
		return 0.54*mean+0.72;
	}
	/***
	 * Returns array of threshold values for input window size away from each arr value. If magnitude is unreasonably 
	 * low(not a step), threshold will be set to an extreme high value to avoid counting peaks in no-step noise areas.
	 * @param arr array of magnitudes to calculate thresholds from
	 * @param windowLength length of "range" of values next to each magnitude value to calculate threshold from
	 * @return array of threshold values for every magnitude value calculated with int windowLength values to the left/right of
	 * magnitude value
	 */
	public static double[] calculateWindow(double[] arr, int windowLength, double minLimit) {
		double[] result=new double[arr.length];
		for (int i = 0; i < arr.length; i++){
			if(i+windowLength < arr.length) {
				if(i-windowLength >= 0) {
					double meanForInterval = calculateMeanInInterval(arr, i-windowLength, i+windowLength);
					double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, i-windowLength, i+windowLength);
					if(magnitudeIsUnreasonblySmall(arr, WINDOW_LENGTH, i, minLimit)){
						result[i]=((meanForInterval+deviationForInterval))*THRESHOLD_MULTIPLE;
					}else{
						result[i] = (meanForInterval+deviationForInterval);
					}
				} else {
					double meanForInterval = calculateMeanInInterval(arr, 0, i+windowLength);
					double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, 0, i+windowLength);
					if(magnitudeIsUnreasonblySmall(arr, WINDOW_LENGTH, i+WINDOW_LENGTH, minLimit)){
						result[i]=((meanForInterval+deviationForInterval))*THRESHOLD_MULTIPLE;
					}else{
						result[i] = (meanForInterval+deviationForInterval);
					}
				}
			} else {
				double meanForInterval = calculateMeanInInterval(arr, i-windowLength, arr.length);
				double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, i-windowLength, arr.length);
				if(magnitudeIsUnreasonblySmall(arr, WINDOW_LENGTH, i-WINDOW_LENGTH, minLimit)){
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
	 * @return true if >80% of values are under the threshold 
	 */
	private static boolean magnitudeIsUnreasonblySmall(double[] arr, int windowLength, int index, double threshold) {
		int underThresholdCounter = 0;
		int totalIterations = 0;
		for(int i = index-windowLength; i < index+windowLength; i++) {
			if(arr[i] < threshold) {
				underThresholdCounter++;
			}
			totalIterations++;
		}
		if((double)underThresholdCounter/(double)totalIterations > 0.8) return true;
		return false;
	}
	
	public static double[] noiseSmoothing(double[] magnitudes, int averageLength) {
		double[] result = new double[magnitudes.length-averageLength];
		result = NoiseSmoothing.generalRunningAverage(magnitudes, averageLength);
		return result;
	}
	
	public static double calculateMagnitude(double x, double y, double z) {
		return Math.sqrt(x*x + y*y + z*z);
	}
	
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
	
	private static double calculateStandardDeviation(double[] arr, double mean) {
		double sum = 0;
		for (int i = 0; i < arr.length; i++) {
			sum += (arr[i] - mean)*(arr[i] - mean);
		}
		return Math.sqrt(sum/(arr.length-1));
	}
	
	private static double calculateStandardDeviationInInterval(double[] arr, double mean, int startInterval, int endInterval) {
		double sum = 0;
		for (int i = startInterval; i < endInterval; i++) {
			sum += (arr[i] - mean)*(arr[i] - mean);
		}
		return Math.sqrt(sum/(double)(endInterval-startInterval));
	}
	
	private static double calculateMean(double[] arr) {
		double sum = 0;
		for(int i = 0; i < arr.length; i++) {
			sum += arr[i];
		}
		return sum/arr.length;
	}
	
	private static double calculateMeanInInterval(double[] arr, int startIndex, int endIndex) {
		double sum = 0;
		for(int i = startIndex; i < endIndex; i++) {
			sum += arr[i];
		}
		return sum/(double)(endIndex-startIndex);
	}
}