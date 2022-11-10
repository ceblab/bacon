from amplify import BinaryPoly, gen_symbols

# イジング変数q_0, q_1を定義
q = gen_symbols(BinaryPoly, 2)

# 目的関数 1 - q_0 * q_1 を定義
f = 1 - q[0] * q[1]

print(f"f = {f}")

from amplify import decode_solution, Solver
from amplify.client import FixstarsClient

# クライアントの設定
client = FixstarsClient()  # Fixstars Optigan
client.parameters.timeout = 1000  # タイムアウト1秒
client.token = "bGezdYx88VisgcbxBgmpP7i0ferFQiZG"  # ローカル環境で使用する場合は、Amplify AEのアクセストークンを入力してください

# ソルバーの構築
solver = Solver(client)  # ソルバーに使用するクライアントを設定

# 問題を入力してマシンを実行
result = solver.solve(f)  # 問題を入力してマシンを実行

for sol in result:  # 複数の解をイテレート

    # sol.values: 決定変数の値（キーをインデックス、値を変数値とする辞書）
    # sol.energy: 目的関数の値（目的関数に決定変数を代入した値）
    solution = decode_solution(q, sol.values)  #  変数配列qをsol.valuesでデコード

    print(f"result: {q} = {solution} (f = {sol.energy})")