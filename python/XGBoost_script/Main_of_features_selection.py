# -- coding:utf-8 --
import pandas as pd
import numpy as np
from utils_features_selection import *
# reload(xgboost)
from xgboost import plot_tree
from sklearn import linear_model
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree
from datetime import datetime
import os
import pydotplus
from save_log import save_print_to_file
from dtreeviz.trees import *
# python matplotlib PDF 不断字
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['figure.figsize'] = 20,20
date_string = datetime.now().strftime("%Y_%m_%d")
#inpath = '/Users/ranpeng/Desktop/胆管癌/data/2020-09-04/XGBoost_452genes.csv'
#df = pd.read_csv(inpath, skiprows=0)
# inpath = '/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/results/XGBoost/XGBoost_select.csv'
# df = pd.read_csv(inpath, skiprows=0)
## features selection part
def features_selection(outpath = None):
    ## 读取375的数据
    X_data_all_features = df.drop('Feature', axis=1)
    Y_data = df['Feature']
    x_col = df.columns[1:]
    # name_dict = {'乳酸脱氢酶':'Lactate dehydrogenase (LDH)','淋巴细胞(%)':'Lymphocytes(%)','超敏C反应蛋白':'High-sensitivity C-reactive protein (hs-CRP)',
    #          '钠':'Sodium','氯':'Chlorine','国际标准化比值':'International Normalized Ratio (INR)','嗜酸细胞(#)':'Eosinophils(#)',
    #          '嗜酸细胞(%)':'Eosinophils(%)','单核细胞(%)':'Monocytes(%)','白蛋白':'Albumin'}
    #
    # 构建一个dataframe用于存储特征的重要程度信息
    import_feature = pd.DataFrame()
    import_feature['col'] = x_col
    import_feature['xgb'] = 0
    # 重复100次试验
    for i in range(100): # 50,150
        #每次试验将375数据随机划分0.7训练集和0.3测试集，注意随机random_state=i
        ## 注明：此方法原因是由于可获得的样本量较少，为了产生不同的训练样本集，使得特征的重要度排序更为稳定，从而选择了这样一种方式。
        ## 通过每次不同的随机种子产生不同的样本，从而达到一定程度上的抑制少量样本的异常对特征的重要度带来的影响。
        x_train, x_test, y_train, y_test = train_test_split(X_data_all_features, Y_data, test_size=0.3, random_state=i)
        #定义模型超参数
        # param1
        # model = xgb.XGBClassifier(
        #         max_depth=3
        #         ,learning_rate=0.2
        #         ,reg_lambda=1
        #         ,n_estimators=200
        #         ,subsample = 0.8
        #         ,colsample_bytree = 0.9)

        # param2
        model = xgb.XGBClassifier(max_depth=3,
                                  learning_rate=0.2,
                                  reg_lambda=1,
                                  n_estimators=150,
                                  subsample=0.9, # original 0.9
                                  colsample_bytree=0.9) # original 0.9
        #模型拟合
        model.fit(x_train, y_train)
        #累加特征重要程度
        import_feature['xgb'] = import_feature['xgb']+model.feature_importances_/100
    # 按照特征重要程度，降序排列
    import_feature = import_feature.sort_values(axis=0, ascending=False, by='xgb')
    print('Top 10 features:')
    print(import_feature.head(10))
    # Sort feature importances from GBC model trained earlier
    # 按照特征重要程度的位置信息
    indices = np.argsort(import_feature['xgb'].values)[::-1]
    #获取前10个重要的特征位置
    Num_f = 12
    indices = indices[:Num_f]

    # Visualise these with a barplot
    # plt.subplots(dpi=400,figsize=(12, 10))
    plt.subplots(figsize=(12, 10))
    # g = sns.barplot(y=list(name_dict.values())[:Num_f], x = import_feature.iloc[:Num_f]['xgb'].values[indices], orient='h') #import_feature.iloc[:Num_f]['col'].values[indices]
    g = sns.barplot(y=import_feature.iloc[:Num_f]['col'].values[indices], x = import_feature.iloc[:Num_f]['xgb'].values[indices], orient='h') #import_feature.iloc[:Num_f]['col'].values[indices]
    g.set_xlabel("Relative importance",fontsize=18)
    g.set_ylabel("Features",fontsize=18)
    g.tick_params(labelsize=14)
    sns.despine()
    plt.savefig(outpath + 'Top 10 feature_importances.pdf')
    # plt.show()
    # g.set_title("The mean feature importance of XGB models");
    # 获取前10重要特征的重要性数值
    import_feature_cols= import_feature['col'].values[:20]
    import_feature.to_csv(outpath + "Top10_features.csv")

    # 画特征金字塔
    num_i = 1
    val_score_old = 0
    val_score_new = 0
    while val_score_new >= val_score_old:
        val_score_old = val_score_new
        # 按重要程度顺序取特种
        x_col = import_feature_cols[:num_i]
        print(x_col)
        X_data = X_data_all_features[x_col]#.values
        ## 交叉验证
        print('5-Fold CV:')
        acc_train, acc_val, acc_train_std, acc_val_std = StratifiedKFold_func_with_features_sel(X_data.values,Y_data.values)
        print("Train AUC-score is %.4f ; Validation AUC-score is %.4f" % (acc_train,acc_val))
        print("Train AUC-score-std is %.4f ; Validation AUC-score-std is %.4f" % (acc_train_std,acc_val_std))
        val_score_new = acc_val
        num_i += 1

    print('Selected features:',x_col[:-1])
    Selected_features = x_col[:-1]
    return list(x_col[:-1])

def single_tree(cols=None):
    cols.insert(0, 'Feature')
    selected_features_data = df[cols]
    selected_features_data.to_csv(outpath + date_string + '_selected_features_data.csv', index=False)
    df_single = selected_features_data
    y_col = df_single['Feature']
    x_col = df_single.drop('Feature', axis = 1)

    # 获取351病人的三特征数据
    x_np = x_col.values
    # 获取351病人的标签数据
    y_np = y_col.values
    # 获取110病人的三特征数据
    #x_test = data_pre_df[x_col].values
    # 在351病人上划分训练集和验证集，此时110视为测试集
    X_train, x_test, y_train, y_test = train_test_split(x_np, y_np, test_size=0.3, random_state=42)
    #限定单树xgb模型
    model = xgb.XGBClassifier(
        max_depth=5,
        n_estimators=1,
    )
    model.fit(X_train,y_train)

    #训练集混淆矩阵
    pred_train = model.predict(X_train)
    #show_confusion_matrix(y_train, pred_train)
    LABELS = ['Good_circulation','Bad_circulation']
    matrix = metrics.confusion_matrix(y_train, pred_train)
    # plt.figure(dpi=400,figsize=(4.5, 3))
    plt.figure(figsize=(4.5, 3))
    sns.heatmap(matrix,
                cmap='coolwarm',
                linecolor='white',
                linewidths=1,
                xticklabels=LABELS,
                yticklabels=LABELS,
                annot=True,
                fmt='d')
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    #plt.show()
    plt.savefig(outpath + date_string + '_train_CM.pdf')
    print(classification_report(y_train, pred_train))
    #plt.savefig('2020-7-20_train.pdf')


    #验证集混淆矩阵
    #pred_val = model.predict(X_val)
    #show_confusion_matrix(y_val, pred_val)
    #print(classification_report(y_val, pred_val))
    #测试集混淆矩阵

    pred_test = model.predict(x_test)
    print('True test label:',y_test)
    print('Predict test label:',pred_test.astype('int32'))
    #show_confusion_matrix(y_test, pred_test)
    LABELS = ['Good_circulation','Bad_circulation']
    matrix = metrics.confusion_matrix(y_test, pred_test)
    # plt.figure(dpi=400,figsize=(4.5, 3))
    plt.figure(figsize=(4.5, 3))
    sns.heatmap(matrix,
                cmap='coolwarm',
                linecolor='white',
                linewidths=1,
                xticklabels=LABELS,
                yticklabels=LABELS,
                annot=True,
                fmt='d')
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    #plt.show()
    plt.savefig(outpath + date_string + '_test_CM.pdf')
    print(classification_report(y_test, pred_test))


    #单树可视化
    features = x_col.columns
    ceate_feature_map(features)
    model_single = xgb.XGBClassifier(
        max_depth=6,
        n_estimators=1,
    )
    model_single.fit(X_train,y_train)
    # model_single.booster().dump_model(outpath + date_string + '_xgb_model.txt')
    graph = xgb.to_graphviz(model_single, fmap='xgb.fmap', num_trees=0, **{'size': str(10)})
    graph.render(filename='single-tree.dot')
    graph.save(outpath + date_string + '_single_tree.dot')

    plt.figure(dpi=300,figsize=(8,6))
    plot_tree(model_single, fmap='xgb.fmap')
    # plot_tree(model)
    # plt.show()
    plt.savefig(outpath + date_string + '_single-tree.pdf')

    # classifier = model  # limit depth of tree
    # # iris = load_iris()
    # # classifier.fit(iris.data, iris.target)

    # viz = dtreeviz(classifier,
    #             x_col,
    #             y_col,
    #             target_name='variety',
    #             feature_names=df.Feature,
    #             class_names=["good", "bad"]  # need class_names for classifier
    #             )

    # viz.view()

    # save split tree graph
    # dot_data = tree.export_graphviz(
    #     model,
    #     out_file=None,
    #     feature_names=df.Feature,
    #     class_names=["good", "bad"],
    #     filled=True,
    #     rounded=True,
    #     special_characters=True,
    # )
    # graph2 = pydotplus.graph_from_dot_data(dot_data)
    # graph2.write_pdf(outpath + date_string + "_tree_split.pdf")


def compare_with_other_method(sub_cols=None, auc_type = 'Train'):

    # 读取351数据集（从375中删除sub_cols全为空的样本得到）
    # x_np,y_np,x_col = data_read_and_split(is_dropna=True,sub_cols=sub_cols)
    # sub_cols.insert(0, 'Feature')
    selected_features_data = df[sub_cols]
    print(selected_features_data)
    # selected_features_data.to_csv(outpath + date_string + '_selected_features_data.csv', index=False)
    compare_model_data = selected_features_data
    x_np = compare_model_data.drop('Feature', axis=1)
    y_np = compare_model_data['Feature']
    x_col = compare_model_data.columns[1:]
    #为了si的图4说明问题。如果是5折，画不出SI 图4
    X_train, X_val, y_train, y_val = train_test_split(x_np, y_np, test_size=0.3, random_state=6)

    #定义全特征下的比对方法
    xgb_n_clf = xgb.XGBClassifier(max_depth=4, # orignal 4
                                  learning_rate=0.2, # orignal 0.3
                                  reg_lambda=1,
                                  n_estimators=150,
                                  subsample=0.9,
                                  colsample_bytree=0.9,
                                  random_state=0)
    tree_clf = tree.DecisionTreeClassifier(random_state=0,max_depth=4) #random_state=0,之前没加
    RF_clf1 = RandomForestClassifier(random_state=0,n_estimators=150,max_depth=4,)
    LR_clf = linear_model.LogisticRegression(random_state=0,C=1,solver='lbfgs')
    LR_reg_clf = linear_model.LogisticRegression(random_state=0,C=0.1, solver='lbfgs')

    fig = plt.figure(dpi=400,figsize=(16, 12))

    Num_iter = 100

    i = 0
    labels_names = []
    Moodel_name = ['Multi-tree XGBoost with all features',
                   'Decision tree with all features',
                   'Random Forest with all features',
                   'Logistic regression with all features with regularization parameter = 1 (by default)',
                   'Logistic regression with all features with regularization parameter = 10',]
    # plot compared multi-model AUC
    AUC_test = ['Train', 'Test']
    for model in [xgb_n_clf,tree_clf,RF_clf1,LR_clf,LR_reg_clf]:
        print('Model:'+Moodel_name[i])
        #以f1的评价方式来k折
        acc_train, acc_val, acc_train_std, acc_val_std = StratifiedKFold_func(x_np.values, y_np.values,Num_iter,model, score_type ='f1')
        #print('F1-score of Train:%.6f with std:%.4f \nF1-score of Validation:%.4f with std:%.6f '%(acc_train,acc_train_std,acc_val,acc_val_std))
        # 以auc的评价方式来k折
        acc_train, acc_val, acc_train_std, acc_val_std = StratifiedKFold_func(x_np.values, y_np.values,Num_iter,model, score_type ='auc')
        print('AUC of Train:%.6f with std:%.4f \nAUC of Validation:%.6f with std:%.4f '%(acc_train,acc_train_std,acc_val,acc_val_std))

        #为了画si图4
        model.fit(X_train,y_train)
        pred_train_probe = model.predict_proba(X_train)[:,1]
        pred_val_probe = model.predict_proba(X_val)[:,1]
        if auc_type == "Test":
            plot_roc(y_val, pred_val_probe,Moodel_name[i],fig,labels_names,i) # 为了画si图4中的test
        elif auc_type == "Train":
            plot_roc(y_train, pred_train_probe,Moodel_name[i],fig,labels_names,i) # 为了画si图4 train
        print('AUC socre:',roc_auc_score(y_val, pred_val_probe))

        i = i+1

    ## 三特征的单树模型对比
    #x_np_sel = x_np[sub_cols] #选择三特征
    ## 划分数据集是为了单树的单次训练并生成AUC图，划分方式和之前保存一致。
    X_train, X_val, y_train, y_val = train_test_split(x_np, y_np, test_size=0.3, random_state=6)

    #为了三特征的模型对比
    xgb_clf = xgb.XGBClassifier(
        max_depth=4,
        n_estimators=1,
        random_state=0,
    )

    tree_clf = tree.DecisionTreeClassifier(random_state=0,max_depth=3)
    RF_clf2 = RandomForestClassifier(random_state=0,n_estimators=1,max_depth=3,)

    #i = 0
    feature_num = len(x_np.columns)
    Moodel_name = ['Single-tree XGBoost with {} features'.format(feature_num),
                'Decision tree with {} features'.format(feature_num),
                'Random Forest with a single {} constraint with three features'.format(feature_num),]
    for model in [xgb_clf,tree_clf,RF_clf2]:
        print('Model'+Moodel_name[i-5])
        #f1的结果
        acc_train, acc_val, acc_train_std, acc_val_std = StratifiedKFold_func(x_np.values, y_np.values,Num_iter,model, score_type ='f1')
        print('F1-score of Train:%.6f with std:%.4f \nF1-score of Validation:%.4f with std:%.6f '%(acc_train,acc_train_std,acc_val,acc_val_std))
        #auc的结果
        acc_train, acc_val, acc_train_std, acc_val_std = StratifiedKFold_func(x_np.values, y_np.values,Num_iter,model, score_type ='auc')
        print('AUC of Train:%.6f with std:%.4f \nAUC of Validation:%.6f with std:%.4f '%(acc_train,acc_train_std,acc_val,acc_val_std))

        model.fit(X_train,y_train)
        pred_train_probe = model.predict_proba(X_train)[:,1]    # 为了画si图4中的train
        pred_val_probe = model.predict_proba(X_val)[:,1]    # 为了画si图4中的test
        if auc_type == "Test":
            plot_roc(y_val, pred_val_probe,Moodel_name[i-5],fig,labels_names,i) # 为了画si图4中的test
        elif auc_type == "Train":
            plot_roc(y_train, pred_train_probe,Moodel_name[i-5],fig,labels_names,i)# 为了画si图4中的train
        print('AUC socre:',roc_auc_score(y_val, pred_val_probe))

        i = i+1

    plt.plot([0,1],[0,1],'r--')
    plt.legend(loc='lower right', fontsize=14)
    if auc_type == "Test":
        plt.savefig(outpath + date_string + '_AUC_test.pdf')
    elif auc_type == "Train":
        plt.savefig(outpath + date_string+ '_AUC_train.pdf')
    # plt.show()

if __name__ == '__main__':

    ## 特征筛选
    inpath = '/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/script/results/XGBoost/Model_blending/Model_blending.csv'
    df = pd.read_csv(inpath, skiprows=0)
    # df.drop(['Lateral_branch_blood_supply'], axis=1, inplace=True)
    outpath = r'/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/script/results/XGBoost/Model_blending/output1_md4_lr0.2/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    # Save the output parameters
    save_print_to_file(outpath)
    # Feature selection
    # selected_cols = features_selection(outpath)
    # selected_cols = ['Basophil_count','Plasma_fibrinogen','Normal_control_PT','Medical_History_Year','Chlorine','Serum_HDL_cholesterol']
    selected_cols = df.columns.tolist()[1:]
    print(f"The select features are: {selected_cols}")
    single_tree(cols=selected_cols)
    # Compare Method
    print('Compare with other methods')
    AUC_test = ['Train', 'Test']
    for auc_type in AUC_test:
        compare_with_other_method(sub_cols=selected_cols, auc_type=auc_type)
