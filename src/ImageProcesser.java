import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import javax.imageio.ImageIO;


public class ImageProcesser {
	
	
	
	
	public int[][]  getVimage(BufferedImage bf){
		int weight = bf.getWidth();
		int height = bf.getHeight();
		int [][]sgray = new int [weight][height];
		 for(int i=0; i<weight; i++) {  		 
	         for(int j=0; j<height; j++) {  
	             int rgb = bf.getRGB(i, j);  
	             int r = (rgb & 0xff0000) >> 16;  
	             int g = (rgb & 0xff00) >> 8;  
	             int b = (rgb & 0xff);     
	             int rgb1 = 0xff000000;
	             rgb1 = rgb1^(r<<16);
	             rgb1 = rgb1^(g<<8);
	             rgb1 = rgb1^(b);
	             int gray = (int)(r * 0.3 + g * 0.59 + b * 0.11);    //计算灰度值  ,其实计算彩色图像的计算公式
	             int mygray = (gray >> 16) | (gray >> 8) | gray;
	            sgray[i][j] = mygray;
	         }  
	     } 
		return sgray;
	}
	
	
	
	
	
	
	
	
	
  public void midFilter(int[][] vimage){
	 // ArrayList arrayList = new ArrayList();
		//BufferedImage grayImage = new BufferedImage(bf.getWidth(),bf.getHeight(), BufferedImage.TYPE_BYTE_GRAY); 
		//int weight = bf.getWidth();
		//int height = bf.getHeight();
		//int [][]sgray = new int [weight][height];
		//int [][]vimage = new int [weight][height];
		int[] temp = new int[9];
		int[][]sgray;
		sgray = vimage;
/*		 for(int i=0; i<weight; i++) {  		 
	         for(int j=0; j<height; j++) {  
	             int rgb = bf.getRGB(i, j);  
	             int r = (rgb & 0xff0000) >> 16;  
	             int g = (rgb & 0xff00) >> 8;  
	             int b = (rgb & 0xff);     
	             int rgb1 = 0xff000000;
	             rgb1 = rgb1^(r<<16);
	             rgb1 = rgb1^(g<<8);
	             rgb1 = rgb1^(b);
	             int gray = (int)(r * 0.3 + g * 0.59 + b * 0.11);    //计算灰度值  ,其实计算彩色图像的计算公式
	             int mygray = (gray >> 16) | (gray >> 8) | gray;
	            sgray[i][j] = mygray;
	         }  
	     } */ 
	     for(int i=0; i<256-2; i++) {               //中值过滤，只处理了出边界以外的像素
	    	 for(int j=0; j<256-2; j++) { 
	        	 temp[0] = sgray[i][j];
	        	 temp[1] = sgray[i+1][j];
	        	 temp[2] = sgray[i+2][j];
	        	 
	        	 
	        	 temp[3] = sgray[i][j+1];
	        	 temp[4] = sgray[i+1][j+1];
	        	 temp[5] = sgray[i+2][j+1];
	        	 
	        	 
	        	 temp[6] = sgray[i][j];
	        	 temp[7] = sgray[i+1][j+2];
	        	 temp[8] = sgray[i+2][j+2];
	        	 Arrays.sort(temp);
	        	vimage[i+1][j+1]=temp[4]; 
	        	 
	         }
	     }
	     
	 //    arrayList.add(vimage); 
	     
	     
/*	     for(int i=0; i<weight-2; i++) {  
	    	 for(int j=0; j<height-2; j++) { 
	    		 
	    		 int r = vimage[i][j];
	    		 int g = vimage[i][j];
	    		 int b = vimage[i][j];
	    		 int rgb1 = 0xff000000;
	    		 rgb1 = rgb1^(r<<16);
	               rgb1 = rgb1^(g<<8);
	               rgb1 = rgb1^(b);
	    		 
	               vimage[i][j] = rgb1;
	    		 
	    		 grayImage.setRGB(i, j, vimage[i][j]);
	    	 }
		 
	     }
	   
	     arrayList.add(grayImage);
	     return arrayList;
	     */
  }
  
  
  
  public int[][] normalize(int [][] vimage){        //统一化
	  int [] histogram = new int[256];
	 // System.out.println(vimage.length);
	  for(int i=1;i<255;i++){
		  for(int j=1;j<255;j++){
			  try{
			//	  System.out.println(vimage[i][j]+""+i+""+j);
			  histogram[vimage[i][j]]++;
			  //System.out.println(vimage[i][j]);
			  }catch(Exception e){
				  
			//	System.out.println("wrong");
			  }
		  }
	  }
	  
	  double dMean = 0;
	  for(int i=1; i<255; i++)
			dMean += i*histogram[i];
		dMean = (dMean/(256*256));
		
		double dSigma = 0;
		for(int i=0; i<255; i++)
			dSigma += histogram[i]*(i-dMean)*(i-dMean);
		dSigma /= (256*256);
		dSigma = Math.sqrt(dSigma);

		double dMean0 = 128, dSigma0 = 128;
		double dCoeff = dSigma0/dSigma;

		for(int i=0; i<256; i++)
		{
			for(int j=0; j<256; j++)
			{
				double dVal = vimage[i][j];
				dVal = dMean0 + dCoeff*( dVal - dMean0);
				if(dVal<0)
					dVal = 0;
				else if (dVal>255)
					dVal = 255;

				vimage[i][j] = (int)dVal;
			}
		}
		
		return vimage;
  }
  
  public void inverseImage(int[][] vimage){
	  
	  for(int i=0; i<256; i++)
		{
			for(int j=0; j<256; j++)
			{
				vimage[i][j]	= (0xFF - vimage[i][j]);
				
			}
		}
  }
  
  
  public static void convertRGB(int[][] vimage){
	    for(int i=0; i<256-2; i++) {  
	    	 for(int j=0; j<256-2; j++) { 
	    		 
	    		 int r = vimage[i][j];
	    		 int g = vimage[i][j];
	    		 int b = vimage[i][j];
	    		 int rgb1 = 0xff000000;
	    		 rgb1 = rgb1^(r<<16);
	               rgb1 = rgb1^(g<<8);
	               rgb1 = rgb1^(b);
	    		 
	               vimage[i][j] = rgb1;

	    	 }
		 
	     }
	   
	  
  }
  public static void setRGB(int [][] vimage,BufferedImage grayImage){
      
      for(int i=0; i<256-2; i++) {  
	    	 for(int j=0; j<256-2; j++) { 
	    		 
	    		 grayImage.setRGB(i, j, vimage[i][j]);
	    		 
	    	 }
      }
  }
  
  
  
  
  public void directionImage(int[][] vimage,double [][] direction){
	int  SEMIBLOCKSIZE = 7;
	//double[][] direction = new double[256][256];
	int [][]dx = new int[SEMIBLOCKSIZE*2+1][SEMIBLOCKSIZE*2+1];
	int [][] dy = new int[SEMIBLOCKSIZE*2+1][SEMIBLOCKSIZE*2+1];
    double	fx, fy;


		for(int y=SEMIBLOCKSIZE+1; y<256-SEMIBLOCKSIZE-1; y++)
		{
			for(int x=SEMIBLOCKSIZE+1; x<256-SEMIBLOCKSIZE-1; x++)
			{
				// calc dx and dy
				for(int j=0; j<SEMIBLOCKSIZE*2+1; j++)
				{
					for(int i=0; i<SEMIBLOCKSIZE*2+1; i++)
					{
						//dx[i][j] = long(ucImg[(y+j-SEMIBLOCKSIZE)*lWidth + x+i-SEMIBLOCKSIZE] - ucImg[(y+j-SEMIBLOCKSIZE)*lWidth + x+i-SEMIBLOCKSIZE-1]);
						//dy[i][j] = long(ucImg[(y+j-SEMIBLOCKSIZE)*lWidth + x+i-SEMIBLOCKSIZE] - ucImg[(y+j-SEMIBLOCKSIZE-1)*lWidth + x+i-SEMIBLOCKSIZE]);
					    dx[i][j] = vimage[y+j-SEMIBLOCKSIZE][x+i-SEMIBLOCKSIZE] - vimage[y+j-SEMIBLOCKSIZE][x+i-SEMIBLOCKSIZE-1];
					    dy[i][j] = vimage[y+j-SEMIBLOCKSIZE][x+i-SEMIBLOCKSIZE] - vimage[y+j-SEMIBLOCKSIZE][x+i-SEMIBLOCKSIZE-1];
					
					}
				}

				// calc direciton
				fx = 0.0;
				fy = 0.0;
				for(int j=0; j<SEMIBLOCKSIZE*2+1; j++)
				{
					for(int i=0; i<SEMIBLOCKSIZE*2+1; i++)
					{
						fx += 2*dx[i][j]*dy[i][j];
						fy += (dx[i][j]*dx[i][j]-dy[i][j]*dy[i][j]);
					}
				}

				// to angle
				direction[y][x] = Math.atan2(fx, fy);
			}
		}

  }
  
  
 public void  dircLowPass(double [][] direction,double [] directionLow ){
        int DIR_FILTER_SIZE = 2;
		int		blocksize = 2*DIR_FILTER_SIZE + 1;
/*	    float* filter = NULL;
	    float* phix   = NULL;
	    float* phiy   = NULL;
	    float* phi2x  = NULL;
	    float* phi2y  = NULL;*/
		int	imgsize = 256*256;
		
	double []	filter = new double[blocksize*blocksize];
	double [] 	phix = new double[imgsize];
	double[]	phiy = new double[imgsize];
	double[]	phi2x = new double[imgsize];
	double[]	phi2y = new double[imgsize];

	

		double	tempSum;
		// filter mask
		tempSum = 0.0;
		for(int y=0; y<blocksize; y++)
		{
			for(int x=0; x<blocksize; x++)
			{
				filter[y*blocksize+x] = (int)(blocksize - (Math.abs(DIR_FILTER_SIZE-x)+Math.abs(DIR_FILTER_SIZE-y)));
				tempSum  += filter[y*blocksize+x];
			}
		}

		for(int y=0; y<blocksize; y++)
		{
			for(int x=0; x<blocksize; x++)
			{
				filter[y*blocksize+x] /= tempSum;
			}
		}

		for(int y=0; y<256; y++)
		{
			for(int x=0; x<256; x++)
			{
				phix[y*256+x] = Math.cos(direction[y][x]);
				phiy[y*256+x] = Math.sin(direction[y][x]);
			}
		}

		// low pass for phi
	//	memset(phi2x, 0, sizeof(float)*imgsize);
		//memset(phi2y, 0, sizeof(float)*imgsize);
		double nx,ny;
		int	val;
	    for (long y = 0; y < 256-blocksize; y++)
		{
			for (long x = 0; x < 256-blocksize; x++)
			{
				nx = 0.0;
				ny = 0.0;
				for (int j = 0; j < blocksize; j++)
				for (int i = 0; i < blocksize; i++)
				{
					val = (int)((x+i)+(j+y)*256);
					nx += filter[j*blocksize+i]*phix[val];
					ny += filter[j*blocksize+i]*phiy[val];
				}
				val = (int)(x+y*256);
				phi2x[val] = nx;
				phi2y[val] = ny;
			}
		}

	    for (long y = 0; y < 256-blocksize; y++)
		{
			for (long x = 0; x < 256-blocksize; x++)
			{
				val = (int)(x+y*256);
	            directionLow[val] = Math.atan2(phi2y[val], phi2x[val])*0.5;
			}
		}
 }
  
  
  
  
  
  
  public static void main(String[] args)throws Exception {
	String path = "F:/kkllor 毕设/新建文件夹/FingerPrintVerify/FingerPrintVerify/2_2.bmp";
	BufferedImage bf = ImageIO.read(new File(path));	 
	  File newFile = new File("F:/kkllor 毕设/新建文件夹/FingerPrintVerify/FingerPrintVerify/1_0_6.bmp"); 
	    if(!newFile.exists()){
	    	newFile.createNewFile();
	    }
	    ImageProcesser ip = new ImageProcesser();
	   // ArrayList arrayList = ip.midFilter(bf);
	    int [][] vimage = ip.getVimage(bf);
	  //  ip.midFilter(vimage);
		BufferedImage grayImage = new BufferedImage(bf.getWidth(),bf.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
		ip.normalize(vimage);
		double[][]direction = new double[256][256];
		double[]direction1 = new double[256*256];
		
		ip.directionImage(vimage, direction);
		ip.dircLowPass(direction, direction1);
		
		
		
		
		ip.inverseImage(vimage);
		
		
		Gabor gabor = new Gabor();
		
		int [][] afterVimage = gabor.ImageEnhance(vimage, direction1);
		
		
       /* ip.convertRGB(vimage);
        ip.setRGB(vimage, grayImage);*/
		 ip.convertRGB(afterVimage);
	        ip.setRGB(afterVimage, grayImage);
	   ImageIO.write(grayImage, "bmp", newFile); 
	   
	   
}
}
