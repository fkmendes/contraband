package testdrivers;

import java.util.Arrays;

import contraband.ManyValuesOneContTrait;

public class ManyValuesOneContTraitTestDriver {

	public static void main(String[] args) {
		String oneTraitValues = "sp1=-0.65187041,2.04994431,0.70321432,0.61298809,0.04211039|sp2=0.9181909,1.7149411,-0.2674609,0.8641253,1.4442747";
		ManyValuesOneContTrait oneTrait = new ManyValuesOneContTrait();
		oneTrait.initByName("traitValues", oneTraitValues);
		
		System.out.println(Arrays.toString(oneTrait.getSample("sp2")));
	}
	
}
