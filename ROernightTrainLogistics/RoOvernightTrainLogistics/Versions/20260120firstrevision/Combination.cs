using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RoOvernightTrainLogistics
{
    class Combination
    {
        public static List<List<T>> Combine<T>(List<T> elements, int k)
        {
            List<List<T>> combinations = new List<List<T>>();
            Combine(elements, k, 0, new List<T>(), combinations);
            return combinations;
        }

        private static void Combine<T>(List<T> elements, int k, int start, List<T> currentCombination, List<List<T>> combinations)
        {
            if (k == 0)
            {
                combinations.Add(new List<T>(currentCombination));
                return;
            }

            for (int i = start; i <= elements.Count - k; i++)
            {
                currentCombination.Add(elements[i]);
                Combine(elements, k - 1, i + 1, currentCombination, combinations);
                currentCombination.RemoveAt(currentCombination.Count - 1);
            }
        }

        //public static void Main(string[] args)
        //{
        //    List<int> numbers = new List<int> { 1, 2, 3, 4 };

        //    int k = 3; // 选择3个数进行组合
        //    List<List<int>> combinations = Combine(numbers, k);

        //    Console.WriteLine("组合（选择" + k + "个数）：");
        //    foreach (List<int> combination in combinations)
        //    {
        //        Console.WriteLine(string.Join(", ", combination));
        //    }
        //}
    }
}
