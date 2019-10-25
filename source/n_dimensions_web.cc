#include "web/Document.h"
#include "web/emfunctions.h"
// #include "web/d3/d3_init.h"
#include "web/d3/selection.h"
#include "web/d3/svg_shapes.h"
#include "web/d3/scales.h"
#include "tools/Graph.h"
#include "n_dimensions.h"

namespace UI = emp::web;

EvoConfig config;

class WebController {
    NDimSim sim;
    UI::Document doc;
    emp::WeightedGraph trans_prob_graph;
   // D3::LinearScale edge_width;

    public:
    WebController(EvoConfig & in_config) : 
        sim(in_config), doc("emp_base") {
        doc << "<svg id='anim_canvas'></svg>";
        int N_GENOTYPES = in_config.N_GENOTYPES();
        trans_prob_graph.Resize(N_GENOTYPES);
        auto mut_rates = sim.GetMutRates();
        double max = 0;
        for (int i = 0; i < N_GENOTYPES; i++) {
            for (int j = 0; j < N_GENOTYPES; j++) {
                double w = mut_rates[i].GetWeight(j);
                if (w) {
                    trans_prob_graph.AddEdge(i, j, w);
                    if (w > max) {
                        max = w;
                    }
                }
            }
        }

        // edge_width.SetDomain(emp::array<double, 3>({0, max}));
        // edge_width.SetRange(emp::array<double, 3>({0, 5}));
    }

    void DrawGraph(double radius = 150) {
        D3::Selection s = D3::Select("#anim_canvas");
        D3::Selection defs = D3::Select("#arrow_defs");

        std::function<std::string(std::string)> MakeArrow = [&defs](std::string color) {
            std::string id = "arrow_" + color;
            emp::remove_chars(id, "#,() ");

            defs.Append("svg:marker")
                .SetAttr("id", id)
                .SetAttr("viewBox", "0 -5 10 10")
                .SetAttr("refX", 5) // This sets how far back it sits, kinda
                .SetAttr("refY", 0)
                .SetAttr("markerWidth", 9)
                .SetAttr("markerHeight", 9)
                .SetAttr("orient", "auto")
                .SetAttr("markerUnits", "userSpaceOnUse")
                .Append("svg:path")
                .SetAttr("d", "M0,-5L10,0L0,5")
                .SetStyle("fill", color);
            
            return "url(#" + id + ")";
        };

        double theta = 0;
        double inc = 2*3.14 / trans_prob_graph.GetNodes().size();
        for (emp::Graph::Node n : trans_prob_graph.GetNodes()) {
            double cx = radius + cos(theta) * radius;
            double cy = radius + sin(theta) * radius;
            theta += inc;
            D3::Selection new_node = s.Append("circle");
            new_node.SetAttr("r", 10).SetAttr("cx", cx).SetAttr("cy", cy);//.BindToolTipMouseover(node_tool_tip);;
            // std::cout << n.GetLabel() << std::endl;
            // EM_ASM_ARGS({js.objects[$0].datum(Pointer_stringify($1));}, new_node.GetID(), n.GetLabel().c_str());
        }

        D3::LineGenerator l;
        auto weights = trans_prob_graph.GetWeights();
        for (int i = 0; i < weights.size(); ++i) {
            for (int j = 0; j < weights[i].size(); ++j) {
                // std::cout << weights[i][j] << std::endl;
                if (weights[i][j]) {
                    double cxi = radius + cos(i*inc) * radius;
                    double cxj = radius + cos(j*inc) * radius;
                    double cyi = radius + sin(i*inc) * radius;
                    double cyj = radius + sin(j*inc) * radius;
                    // std::string color = edge_color.ApplyScaleString(weights[i][j]);
                    // std::cout << color << std::endl;
                    emp::array<emp::array<double, 2>, 2> data({emp::array<double, 2>({cxi, cyi}), emp::array<double, 2>({cxj, cyj})});
                    D3::Selection n = s.Append("path");
                    n.SetAttr("d", l.Generate(data))
                    .SetStyle("stroke", "black")
                    // .SetStyle("stroke-width", edge_width.ApplyScale(weights[i][j]))
                    .SetAttr("marker-end", MakeArrow("black"));
                    // .BindToolTipMouseover(edge_tool_tip);

                    // EM_ASM_ARGS({js.objects[$0].datum($1);}, n.GetID(), weights[i][j]);
                    //  std::cout << MakeArrow(color) << std::endl;
                    // s.Append("path").SetAttr("d", "M"+emp::to_string(cxi)+","+emp::to_string(cyi)+"L"+emp::to_string(cxj)+","+emp::to_string(cyj)).SetStyle("stroke", color).SetStyle("stroke-width", weights[i][j]*10);
                }
            }
        }

        // s.SetupToolTip(edge_tool_tip);
        // s.SetupToolTip(node_tool_tip);

    }

};

WebController controller(config);

int main()
{
    emp::web::OnDocumentReady([](){
        controller.DrawGraph();
    });

}
