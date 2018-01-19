package Utils;

import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;

public class TimeFormat {
	public static Date parse(String str){
		DateFormat df = null;
		if(str.contains("-"))
			df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS", Locale.US);
		else if(str.contains("/"))
			df = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		else if(str.length()==26||str.length()==25)
			df = new SimpleDateFormat("MMM dd yyyy hh:mm:ss:SSSa",Locale.US);
		else if(str.length()==17)
			df = new SimpleDateFormat("yyyyMMdd HH:mm:ss");
		else if(str.length()==12)
			df = new SimpleDateFormat("yyyyMMddHHmm");
		try {
			return df.parse(str);
		} catch (ParseException e) {
			e.printStackTrace();
		}
		return null;
	}

	public static void main(String[] args) {
		//Long t = tf.parse("Oct 29 2014 12:00:00:056AM").getTime();
		System.out.println(TimeFormat.parse("201412291731").getTime());
	}
}
