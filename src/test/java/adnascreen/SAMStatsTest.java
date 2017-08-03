package adnascreen;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.json.JSONObject;

import org.junit.Test;

/**
 * Test code for counting read targets after alignment
 * @author mmah
 *
 */
public class SAMStatsTest {
	String singleTargets = "{"
			+ "'1':'1',"
			+ "'2':'2',"
			+ "'3':'3',"
			+ "'4':'4',"
			+ "'5':'5',"
			+ "'6':'6',"
			+ "'7':'7',"
			+ "'8':'8',"
			+ "'9':'9',"
			+ "'10':'10',"
			+ "'11':'11',"
			+ "'12':'12',"
			+ "'13':'13',"
			+ "'14':'14',"
			+ "'15':'15',"
			+ "'16':'16',"
			+ "'17':'17',"
			+ "'18':'18',"
			+ "'19':'19',"
			+ "'20':'20',"
			+ "'21':'21',"
			+ "'22':'22',"
			+ "'X':'X',"
			+ "'Y':'Y',"
			+ "'MT':'MT',"
			+ "}";
	
	@Test
	public void singleTargetTest(){
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("target-test.sam").getPath();
		try{
			IndexAndBarcodeKey key = new IndexAndBarcodeKey("47_2_Q20_Q41");
			JSONObject targetJSON = new JSONObject(singleTargets);
			SAMStats stats = new SAMStats(key, filename, targetJSON, 0);
			
			SampleCounter counter = stats.counter.get(key);
			// this file is constructed to have increasing number of reads with chromosome
			for (int i = 1; i <= 22; i++){
				assertEquals(i, counter.get(String.valueOf(i)));
			}
			assertEquals(1, counter.get("X"));
			assertEquals(2, counter.get("Y"));
			assertEquals(3, counter.get("MT"));
		}
		catch(Exception e){
			fail();
		}
	}
	
	@Test
	public void groupedTargetTest(){
		String groupedTargets = "{'autosome':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],'X':'X','Y':'Y','MT':'MT'}";
		
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("target-test.sam").getPath();
		try{
			IndexAndBarcodeKey key = new IndexAndBarcodeKey("47_2_Q20_Q41");
			JSONObject targetJSON = new JSONObject(groupedTargets);
			SAMStats stats = new SAMStats(key, filename, targetJSON, 0);
			
			SampleCounter counter = stats.counter.get(key);
			assertEquals(253, counter.get("autosome"));
			assertEquals(1, counter.get("X"));
			assertEquals(2, counter.get("Y"));
			assertEquals(3, counter.get("MT"));
		}
		catch(Exception e){
			fail();
		}
	}
	
	@Test
	public void humanTargetTest(){
		String groupedTargets = "{'human':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']}";
		
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("target-test.sam").getPath();
		try{
			IndexAndBarcodeKey key = new IndexAndBarcodeKey("47_2_Q20_Q41");
			JSONObject targetJSON = new JSONObject(groupedTargets);
			SAMStats stats = new SAMStats(key, filename, targetJSON, 0);
			
			SampleCounter counter = stats.counter.get(key);
			assertEquals(259, counter.get("human"));
		}
		catch(Exception e){
			fail();
		}
	}
	
	@Test
	public void combinationTargetTest(){
		String groupedTargets = "{'autosome':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],'X':'X','Y':'Y','MT':'MT',"
				+ "'human':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']}";
		
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("target-test.sam").getPath();
		try{
			IndexAndBarcodeKey key = new IndexAndBarcodeKey("47_2_Q20_Q41");
			JSONObject targetJSON = new JSONObject(groupedTargets);
			SAMStats stats = new SAMStats(key, filename, targetJSON, 0);
			
			SampleCounter counter = stats.counter.get(key);
			assertEquals(259, counter.get("human"));
			assertEquals(253, counter.get("autosome"));
			assertEquals(1, counter.get("X"));
			assertEquals(2, counter.get("Y"));
			assertEquals(3, counter.get("MT"));
		}
		catch(Exception e){
			fail();
		}
	}
	
	@Test
	public void qualityTest(){
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("target-test.sam").getPath();
		try{
			IndexAndBarcodeKey key = new IndexAndBarcodeKey("47_2_Q20_Q41");
			JSONObject targetJSON = new JSONObject(singleTargets);
			SAMStats stats = new SAMStats(key, filename, targetJSON, 30);
			
			SampleCounter counter = stats.counter.get(key);
			assertEquals(0, counter.get("1"));
			assertEquals(1, counter.get("2"));
			assertEquals(1, counter.get("3"));
			assertEquals(2, counter.get("4"));
			assertEquals(1, counter.get("5"));
			assertEquals(4, counter.get("6"));
			assertEquals(1, counter.get("7"));
			assertEquals(2, counter.get("8"));
			assertEquals(2, counter.get("9"));
			assertEquals(6, counter.get("10"));
			assertEquals(1, counter.get("11"));
			assertEquals(1, counter.get("12"));
			assertEquals(6, counter.get("13"));
			assertEquals(4, counter.get("14"));
			assertEquals(8, counter.get("15"));
			assertEquals(7, counter.get("16"));
			assertEquals(5, counter.get("17"));
			assertEquals(9, counter.get("18"));
			assertEquals(6, counter.get("19"));
			assertEquals(5, counter.get("20"));
			assertEquals(10, counter.get("21"));
			assertEquals(2, counter.get("22"));
			assertEquals(1, counter.get("X"));
			assertEquals(0, counter.get("Y"));
			assertEquals(0, counter.get("MT"));
		}
		catch(Exception e){
			fail();
		}
	}
}
