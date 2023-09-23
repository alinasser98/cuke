#include <torch/extension.h>

torch::Tensor sub_add_a_b_c(torch::Tensor obj_a, torch::Tensor obj_b, torch::Tensor obj_c)
{
    auto a = obj_a.accessor<float, 2>();
auto b = obj_b.accessor<float, 2>();
torch::Tensor obj_arr2 = torch::empty({10,10}, at::kFloat);
auto arr2 = obj_arr2.accessor<float, 2>();
for (int _l0 = 0; _l0 < 10; _l0 += 1) {
for (int _l1 = 0; _l1 < 10; _l1 += 1) {
arr2[_l0][_l1] = a[_l0][_l1] + b[_l0][_l1];
} 
} 
auto c = obj_c.accessor<float, 2>();
torch::Tensor obj_arr6 = torch::empty({10,10}, at::kFloat);
auto arr6 = obj_arr6.accessor<float, 2>();
for (int _l2 = 0; _l2 < 10; _l2 += 1) {
for (int _l3 = 0; _l3 < 10; _l3 += 1) {
arr6[_l2][_l3] = arr2[_l2][_l3] - c[_l2][_l3];
} 
} 
return obj_arr6;

}

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
    m.def("run", &sub_add_a_b_c);
}