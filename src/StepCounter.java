import java.util.Arrays;

public class StepCounter {
	private static final int WINDOW_LENGTH = 3;
	private static final int THRESHOLD_LIMIT = 1;
	private static final int MINIMUM_LIMIT = 2;
	private static final int THRESHOLD_MULTIPLE=5;
	public static void main(String[] args) {
		String[] columnNames={"time", "accel-x", "accel-y", "accel-z"};
		CSVData test = new CSVData("data/50StepWalkMale.csv", columnNames, 1);
		test.correctTime(test);
		System.out.println(countSteps(test.getColumn(1),test.getRows(1,test.getNumRows()-1), WINDOW_LENGTH));
	}
	
	private static int countSteps(double[] times, double[][] sensorData, int windowLength) {
		int stepCount = 0;
		double[] arr = new double[times.length];
		arr = calculateMagnitudes(sensorData);
		//arr = NoiseSmoothing.generalRunningAverage(arr, 3);
		double[] thresholds = calculateWindow(arr, WINDOW_LENGTH);
		for(int i = 1; i < arr.length; i++) {
			if (arr[i] > arr[i-1] && arr[i] > arr[i+1]) {
				if(arr[i] > thresholds[i]) {
					stepCount++;
					System.out.println(stepCount + " " + times[i]);
				}
			}
		}
		return stepCount;	
	}	
	
	/***
	 * Returns array with threshold values
	 * @param arr
	 * @param windowLength
	 * @return
	 */
	public static double[] calculateWindow(double[] arr, int windowLength) {
		double[] result=new double[arr.length];
		for (int i = 0; i < arr.length; i++){
			if(i+windowLength < arr.length) {
				if(i-windowLength >= 0) {
					double meanForInterval = calculateMeanInInterval(arr, i-windowLength, i+windowLength);
					double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, i-windowLength, i+windowLength);
					if(belowThreshold(arr, WINDOW_LENGTH, i, MINIMUM_LIMIT)){
						result[i]=((meanForInterval+deviationForInterval))*THRESHOLD_MULTIPLE;
					}else{
						result[i] = (meanForInterval+deviationForInterval);
					}
				} else {
					double meanForInterval = calculateMeanInInterval(arr, 0, i+windowLength);
					double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, 0, i+windowLength);
					if(belowThreshold(arr, WINDOW_LENGTH, i+WINDOW_LENGTH, MINIMUM_LIMIT)){
						result[i]=((meanForInterval+deviationForInterval))*THRESHOLD_MULTIPLE;
					}else{
						result[i] = (meanForInterval+deviationForInterval);
					}
				}
			} else {
				double meanForInterval = calculateMeanInInterval(arr, i-windowLength, arr.length);
				double deviationForInterval = calculateStandardDeviationInInterval(arr, meanForInterval, i-windowLength, arr.length);
				if(belowThreshold(arr, WINDOW_LENGTH, i-WINDOW_LENGTH, MINIMUM_LIMIT)){
					result[i]=((meanForInterval+deviationForInterval))*THRESHOLD_MULTIPLE;
				}else{
					result[i] = (meanForInterval+deviationForInterval);
				}
			}
		}
		return result;
	}
	
	private static boolean belowThreshold(double[] arr, int windowLength, int index, int threshold) {
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
	