// Simplex.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include"Simplex_Model.h"
int main()
{
    //std::cout << "Hello World!\n";
    //Simplex_Reader sr("zhengshu1.txt");
    //Simplex_Model sm;
	//try
	//{
	//	sr.readin_info(sm);
	//	sm.solve();
	//}
	//catch (MyException& e)
	//{
	//	cout << e.what() << endl;
	//}
	
	Integer_Model im("check_input.txt");
	try
	{
		im.solve();
		im.print_res();
	}
	catch (MyException& e)
	{
		cout << e.what() << endl;
	}

	//vector<double> v = { 1.5,2.5,5,6,7.6 };
	//cout << sm.find_vector_second_max_idx(v);
	//_base_matrix matrix(3, 5);
	//matrix.set(0, 0, 1);matrix.set(0, 1, 2);matrix.set(0, 2, 3); matrix.set(0, 3, 4); matrix.set(0, 4, 5);
	//matrix.set(1, 0, 2); matrix.set(1, 1, 2); matrix.set(1, 2, 4); matrix.set(1, 3, 4); matrix.set(1, 4, 4);
	//matrix.set(2, 0, 5); matrix.set(2, 1, 5); matrix.set(2, 2, 5); matrix.set(2, 3, 5); matrix.set(2, 4, 5);
	//matrix.print_form();
	//vector<int> v = { 0,4 };
	//matrix._del_cols(v);
	//cout << "after" << endl;
	//matrix.print_form();

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
